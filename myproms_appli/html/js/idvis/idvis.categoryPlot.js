/*
################################################################################
# idvis.categoryPlot.js    1.2.0                                               #
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
                  CategoryPlot object
******************************************************************/
idvis.categoryPlot = function(cpData) {
	var CP=this;
	CP.chartType='categoryPlot';
	CP.chartID=idvis.lib.registerChart(CP);
	
	/*****************************************************************
	  CategoryElement object (must be declared here because called at categoryPlot instantiation)
	******************************************************************/
	var categoryElement = function(dataset,catIdx,elData,baseValue) {
		//this.dsIndex=dsIdx;
		this.dataset=dataset;
		this.catIndex=catIdx;
		this.value=elData.value*1;
		this.label=elData.label || null;
		this.externalID=elData.id || this.label;
		this.color=(elData.color)? elData.color : (dataset.params.chart.categories[catIdx].color)? dataset.params.chart.categories[catIdx].color : dataset.params.color; // overwritten if dataset.params.elColorRule
		if (dataset.params.properties) {
			this.propValues=[];
			for (let i=0; i<dataset.params.properties.length; i++) {
				this.propValues[i]=(elData.propValues)? elData.propValues[i] : null;
			}
		}
		this.metaInfo=elData.meta || null; // Array of extra text info
		this.baseValue=baseValue || 0;
		this.svg=null;
	};
	categoryElement.prototype = { // used for point only!!!! :-(
		info: function(type) {
			let infoStrg;
			if (this.dataset.params.type==='point') { // forced to 'max' at hover ~line 1272
				infoStrg=this.label || this.dataset.params.chart.categories[this.catIndex].label || 'No label';
				let valueLabel=this.dataset.params.chart.axesTitles[this.dataset.params.axis].shortTitle || 'Value';
				if (type === 'min') {
					infoStrg+='\n'+valueLabel+': '+idvis.formatNumber(this.value); //.toFixed(3);
				}
				else {
					const dataTypeLabel= this.dataset.params.chart.datasetsLabel || 'Series';
					if (this.dataset.params.chart.datasets.length > 1) infoStrg+='\n'+dataTypeLabel+': '+this.dataset.params.name;
					infoStrg+='\n'+valueLabel+': '+idvis.formatNumber(this.value); //.toFixed(3);
				}
			}
			if (this.propValues) {
				for (let i=0; i<this.dataset.params.properties.length; i++) {
					if (this.propValues[i] !== null) infoStrg+='\n'+this.dataset.params.properties[i]+': '+idvis.formatNumber(this.propValues[i]); // no effect if string
				}
			}
			if (this.metaInfo) {
				for (let i=0; i<this.metaInfo.length; i++) {
					infoStrg+='\n'+this.metaInfo[i];
				}
			}
			return infoStrg;
		}
	};
	
	
	/* Configuration variables */
	this.divID=cpData.div;
	this.div=document.getElementById(this.divID);
    var mainDivID=cpData.div;
	const canvasDivID=this.divID+'_canvas',
	      formDivID=this.divID+'_form'; // form always displayed
	this.formDiv=null; // defined at draw()
	this.getDivID=function(){return canvasDivID;};
	this.horizontalCP=(cpData.orientation && cpData.orientation.match('^h'))? true : false; // also used by idvis (default is vertical)
	var convertAxis=function(oldAxis) {
		var axis=(!oldAxis)? 'N1' : oldAxis+''; // N: numerical axes 1/2, C: category axis
		if (axis.match('^C')) {axis=(CP.horizontalCP)? 'Y' : 'X';} // value is an index of a category
		else {
			if (CP.horizontalCP) {axis=(axis.match('2'))? 'X2' : 'X';}
			else {axis=(axis.match('2'))? 'Y2' : 'Y';}
		}
		return axis;
	};
	this.exportAsImage=cpData.exportAsImage; // null or array [button name, image file name, action script path]

//this.convertValue=function(axis,value) {return value;}; // for threshold line only!!!
	var startAt0=(cpData.startAtzero)? true : false; // ignored if values < 0
	var minCatPlotValue={},maxCatPlotValue={};
	this.plotArea=null;
	this.axesTitles={};
	var axesUsed={}, numberOfAxes=0; // Numerical axes only
	this.features=[];
	this.editableThresholds=[]; // list of editable features
	this.noColorEdition=(cpData.noColorEdition===false)? false : true; // true by default

	var categoryOnClick=cpData.categoryOnClick,
		boxOnClick=cpData.boxOnClick,
		barOnClick=cpData.barOnClick;
	//var excludeOnClick=(cpData.excludeOnClick)? cpData.excludeOnClick : null;
	//var outlierOnClick=(cpData.outlierOnClick)? cpData.outlierOnClick : function(){};
	this.pointOnClick=cpData.pointOnClick; // "this" because used by idvis
	var boxPointExclusion=cpData.boxPointExclusion || {}; // {label,color} for boxplot excluded points
	this.showScales=(cpData.showScales===undefined || cpData.showScales===null || cpData.showScales)? true : false;
	this.chartSettings={current:{}}; // Needed by idvis in case editable thresholds
	/******** Chart variables & objects ********/
	//var colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'];
	//var colorIdx=0;
	const scaleGradients={ // down:1->0, up:0->up
		'RB':function(down,up){return 'rgb('+Math.round(100*down)+'%,0%,'+Math.round(100*up)+'%)'}, // Red->Blue
		'RG':function(down,up){return 'rgb('+Math.round(100*down)+'%,'+Math.round(90*up)+'%,0%)'}, // Red->Green
		'Rb':function(down,up){return 'rgb('+Math.round(100*down)+'%,0%,0%)'}, // Red->black

		'GR':function(down,up){return 'rgb('+Math.round(90*up)+'%,'+Math.round(100*down)+'%,0%)'}, // Green->Red
		'GB':function(down,up){return 'rgb(0%,'+Math.round(90*down)+'%,'+Math.round(100*up)+'%)'}, // Green->Blue
		'Gb':function(down,up){return 'rgb(0%,'+Math.round(90*down)+'%,0%)'}, // Green->black

		'BR':function(down,up){return 'rgb('+Math.round(100*up)+'%,0%,'+Math.round(100*down)+'%)'}, // Blue->Red
		'BG':function(down,up){return 'rgb(0%,'+Math.round(90*up)+'%,'+Math.round(100*down)+'%)'}, // Blue->Green
		'Bb':function(down,up){return 'rgb(0%,0%,'+Math.round(100*down)+'%)'}, // Blue->black
		
		'YR':function(down,up){return 'rgb('+Math.round(Math.max(90*down,100*up))+'%,'+Math.round(90*down)+'%,0%)'}, // Yellow->Red
		'YG':function(down,up){return 'rgb('+Math.round(90*down)+'%,'+Math.round(Math.max(90*down,90*up))+'%,0%)'}, // Yellow->Green
		'YB':function(down,up){return 'rgb('+Math.round(90*down)+'%,'+Math.round(90*down)+'%,'+Math.round(100*up)+'%)'}, // Yellow->Blue
		'Yb':function(down,up){return 'rgb('+Math.round(90*down)+'%,'+Math.round(90*down)+'%,0%)'}, // Yellow->black
		
		'bR':function(down,up){return 'rgb('+Math.round(100*up)+'%,0%,0%)'}, // black->Red
		'bG':function(down,up){return 'rgb(0%,'+Math.round(90*up)+'%,0%)'}, // black->Green
		'bB':function(down,up){return 'rgb(0%,0%,'+Math.round(100*up)+'%)'}, // black->Blue
		'bY':function(down,up){return 'rgb('+Math.round(90*up)+'%,'+Math.round(90*up)+'%,0%)'} // black->Yellow
	};
	const defaultGradient='RB'; //'RB';
	
	/*** Categories ***/
	var mainLayout=(cpData.layout && cpData.layout.match('^s'))? 'stacked' : (cpData.layout && cpData.layout.match('^o'))? 'overlay' : 'normal';
	if (mainLayout !== 'normal') {cpData.groups=[[[],mainLayout]];} // overwrites any group info => 1 will be created with all datasets
	var maxLabel='',baseValues=[],binCount=[];
	var histoClass=false,spaceBeforeBins=0,spaceAfterBins=0; // Changed only for histClass when a category feature is outsite distribution values range
	var histoClassOptRange,minHistClassBinValue,maxHistClassBinValue,binSize; // only for histogram & area
	const categoryLabelColor=cpData.categoryLabelColor || '#000';
	this.categories=[];
	if (cpData.categories.list) { // null for histogram
		let maxLabelLength=0;
		for (let i=0; i<cpData.categories.list.length; i++) {
			let cat={
				label:cpData.categories.list[i][0],
				externalID:cpData.categories.list[i][1] || i, //cpData.categories.list[i][0],
				color:cpData.categories.list[i][2], // overwrites dataset color; null usually
				labelColor:cpData.categories.list[i][2] || categoryLabelColor,
				index:i
			};
			cat.popupText=cpData.categories.list[i][3] || cat.label;
			this.categories.push(cat);
			if (maxLabelLength < cat.label.length) {
				maxLabelLength=cat.label.length;
				maxLabel=cat.label;
			}
			baseValues[i]=[]; // needed if stacking
		}
	}
	else { // histogram(s) or area
		histoClass=true;
		startAt0=true;
		this.categories.length=(cpData.categories.bins && cpData.categories.bins.number)? cpData.categories.bins.number : 20;
		for (let i=0; i<this.categories.length; i++) {
			this.categories[i]={index:i};
			baseValues[i]=[]; // needed if stacking
			binCount[i]=0;
		}
		this.histLimits=cpData.histLimits || [null,null];
	}
	
	/*** Scales ***/
	this.scales=cpData.scales; //{map:<property|attribute to use value from>,range:[beg,end](opt)}
	var scalesFocus={};
	var scalesPointShapes={}; // records point shape of 1st dataset matching a size mapping property
	if (this.scales) {
		/* Color */
		if (this.scales.color) {
			scalesFocus.color=[];
			this.scales.color.forEach((rule)=>{
				scalesFocus.color.push(rule.map);
				if (!rule.gradient) rule.gradient=defaultGradient; // default gradient Red->Blue
				if (rule.transform) {
					if (!rule.transform.label) rule.transform.label='transformed('+rule.map+')';
					if (!rule.transform.func) {
						alert('ERROR: No data transformation function provided for "'+rule.map+'" color scale.\nRaw values will be used.');
						rule.transform=null; 
					}
				}
			});
		}
		/* Point size */
		if (this.scales.pointSize) {
			scalesFocus.pointSize=[];
			this.scales.pointSize.forEach((rule)=>{
				scalesFocus.pointSize.push(rule.map);
				if (rule.transform) {
					if (!rule.transform.label) rule.transform.label='transformed('+rule.map+')';
					if (!rule.transform.func) {
						alert('WARNING: No data transformation function provided for "'+rule.map+'" point size scale.\nRaw values will be used.');
						rule.transform=null; 
					}
				}
			});
		}
	}

	/*** Datasets ***/
    this.datasets=[];
	this.datasetsLabel=cpData.datasetsLabel || cpData.dataSetsLabel; // can be undef
	this.multiPointDatasetRule=cpData.multiDatasetRule || cpData.multiPointDatasetRule; // string, updated by idvis if undefined
	this.hasBoxes=0;
	this.hasPoints=0;
	for (let i=0; i<cpData.datasets.length; i++) {
		let set=cpData.datasets[i];
		//if (set.type==null) {alert('!!!!ERROR: No type defined for dataset '+i); return;}
		this.datasets[i]={};
		this.datasets[i].params={};
		this.datasets[i].params.visible=2;
		this.datasets[i].params.chart=CP;
		this.datasets[i].params.index=i;
		this.datasets[i].params.type=set.type || 'bar'; // default
		this.datasets[i].params.name=set.name || 'dataset '+(i+1);
		this.datasets[i].params.externalID=set.id || this.datasets[i].params.index;
		this.datasets[i].params.properties=set.properties || null; // [array of properties]
		/** Color (&Point shape) **/
		if (this.datasets[i].params.type==='point') {
			[this.datasets[i].params.pointShape,this.datasets[i].params.color]=idvis.lib.generateDatasetAES(this,i,set);
		}
		else {
			this.datasets[i].params.color=set.color || idvis.colorList[i % idvis.colorList.length];
		}
		if (set.appliedScales && scalesFocus.color) { // Find matching scale => overwrites elements color defined by this.datasets[i].params.color
			/* Check chart scales */
			let match=false;
			for (let s=0; s<set.appliedScales.length; s++) {
				const scIdx=scalesFocus.color.indexOf(set.appliedScales[s]);
				if (scIdx >= 0) {
					this.datasets[i].params.colorRule=this.scales.color[scIdx];
					match=true;
					break;
				}
			}
			if (match) {
				/* Check dataset properties */
				this.datasets[i].params.colorRule.mapPropIndex=(this.datasets[i].params.properties)? this.datasets[i].params.properties.indexOf(this.datasets[i].params.colorRule.map) : -1;
				if (!this.datasets[i].params.colorRule.mapPropIndex===-1) {
					alert('idvis chart error: Color mapped property "'+set.colorScale+'" not found for dataset "'+this.datasets[i].params.name+'"');
				}
				this.datasets[i].params.colorRule.userRange=(this.datasets[i].params.colorRule.range)? true : false;
			}
		}
		this.datasets[i].params.opacity=set.opacity || 0.75;
		this.datasets[i].params.selectable=false;
		if (this.datasets[i].params.type==='box') {this.hasBoxes++;}
		else if (this.datasets[i].params.type==='point') {
			this.datasets[i].params.selectable=true;
			/** Point size **/
			this.datasets[i].params.pointSize=set.pointSize || 7; // number or {map:<property|attribute to use value from>,range:[beg,end](opt)}
			if (set.appliedScales && scalesFocus.pointSize) { // Find matching scale => overwrites this.datasets[i].params.pointSize
				/* Check chart scales */
				let match=false;
				for (let s=0; s<set.appliedScales.length; s++) {
					const scIdx=scalesFocus.pointSize.indexOf(set.appliedScales[s]);
					if (scIdx >= 0) {
						this.datasets[i].params.pointSizeRule=this.scales.pointSize[scIdx];
						match=true;
						break;
					}
				}
				if (match) {
					/* Check dataset properties */
					this.datasets[i].params.pointSizeRule.mapPropIndex=(this.datasets[i].params.properties)? this.datasets[i].params.properties.indexOf(this.datasets[i].params.pointSizeRule.map) : -1;
					if (this.datasets[i].params.pointSizeRule.mapPropIndex===-1) {
						alert('idvis chart error: Point size mapped property "'+set.pointSizeScale+'" not found for dataset "'+this.datasets[i].params.name+'"');
					}
					else {
						const mappedProperty=this.datasets[i].params.properties[this.datasets[i].params.pointSizeRule.mapPropIndex];
						if (!scalesPointShapes[mappedProperty]) {scalesPointShapes[mappedProperty]=this.datasets[i].params.pointShape;} // [property]=dataset pointShape
					}
					this.datasets[i].params.pointSizeRule.userRange=(this.datasets[i].params.pointSizeRule.range)? true : false;
				}
			}
			this.hasPoints++;
		}
		if (set.connection) {
			let connect={
				//isConnection:true,
				lineType:(set.connection[0].match('^c'))? 'curve' : 'line',
				color:set.connection[1] || this.datasets[i].params.color,
				width:set.connection[2] || 3,
				pattern:set.connection[3] || null, // can be undefined
				opacity:set.connection[4] || 0.7,
				ignoreMissing:set.connection[5] || false,
				svg:null // SVG path
			};
			//this.datasets[i].params.line=connect;
			this.datasets[i].dataConnection=connect;
		}
		this.datasets[i].params.axis=convertAxis(set.axis);
		//this.datasets[i].params.onclick=(set.onclick)? set.onclick : null;
		//this.datasets[i].params.visible=true;
		this.datasets[i].data=[];
		if (this.datasets[i].params.type=='bar') {startAt0=true;}
		if (mainLayout != 'normal') {cpData.groups[0][0].push(this.datasets[i].params.externalID);}
	}
	
	/*** Groups ***/
	var datasetGroups=[],groupFinalIndex=[];
	var gIdx=-1;
	for (let i=0; i<CP.datasets.length; i++) {
		var dsMatched=false;
		if (cpData.groups) {
			GR:for (let g=0; g<cpData.groups.length; g++) {
				for (let s=0;s<cpData.groups[g][0].length; s++) {
					if (cpData.groups[g][0][s]==CP.datasets[i].params.externalID) { // current dataset is matched!
						if (groupFinalIndex[g]==null) { // first dataset of group id matched
							gIdx=datasetGroups.length;
							groupFinalIndex[g]=gIdx;
							datasetGroups[gIdx]={content:[]};
							datasetGroups[gIdx].layout=(!cpData.groups[g][1])? 'normal' : (cpData.groups[g][1].match('^o'))? 'overlay': (cpData.groups[g][1].match('^s'))? 'stacked' : 'normal';
							if (datasetGroups[gIdx].layout=='stacked') {
								datasetGroups[gIdx].axis=CP.datasets[i].params.axis; // only defined for 'stacked'
							}
						}
						else {
							gIdx=groupFinalIndex[g];
							if (datasetGroups[gIdx].axis && datasetGroups[gIdx].axis != CP.datasets[i].params.axis) { // axis conflict for stacked layout!
								//if (!datasetGroups[gIdx].axisAlert) { // alert only once
								//	datasetGroups[gIdx].axisAlert=true;
								//	var axisNumber=(datasetGroups[gIdx].axis.match('1'))? '1' : '2';
								//	alert('WARNING: Axis conflict detected in group #'+gIdx+'. Using axis '+axisNumber+' for all datasets in group.');
								//}
								//CP.datasets[i].params.axis=datasetGroups[gIdx].axis; // overwrite dataset axis with group axis
								alert('ERROR: Axis conflict detected for dataset group #'+(gIdx+1)+'. Drawing has failed!');
								return;
							}
						}
						datasetGroups[gIdx].content.push(i);		
						dsMatched=true;
						break GR;
					}
				}
			}
		}
		if (dsMatched===false) { // dset not grouped
			gIdx=datasetGroups.length;
			datasetGroups[gIdx]={content:[i],axis:CP.datasets[i].params.axis,layout:mainLayout};
		}
		axesUsed[CP.datasets[i].params.axis]=1;
		CP.datasets[i].params.group=gIdx;
	}
	
	/*** Bars & points values ***/
	if (histoClass===false) {
		for (let i=0; i<cpData.datasets.length; i++) {
			let set=cpData.datasets[i];
			if (!set.values) {continue;}
			if (set.type.match('^bar|point')) {
				//let axis=(set.axis && (set.axis+'').match('2'))? axisRoot+'2' : axisRoot; //+'1';
				let axis=convertAxis(set.axis);
				let gIdx=CP.datasets[i].params.group;
				for (let catIdx=0; catIdx<set.values.length; catIdx++) {
					if (!set.values[catIdx]) continue;
					if (baseValues[catIdx][gIdx]==null) {baseValues[catIdx][gIdx]=0;}
					//var elData=data[catIdx].split(':');
					var newEl=new categoryElement(CP.datasets[i],catIdx,set.values[catIdx],baseValues[catIdx][gIdx]);
					CP.datasets[i].data[catIdx]=newEl;

					/* Color range */
					if (CP.datasets[i].params.colorRule && !CP.datasets[i].params.colorRule.userRange) { // => build range based on actual values 
						let rule=CP.datasets[i].params.colorRule,
							mappedValue=newEl.propValues[rule.mapPropIndex];
						if (idvis.isNumber(mappedValue)) {
							if (rule.range) {
								rule.range[0]=(idvis.isNumber(rule.range[0]))? Math.min(rule.range[0],mappedValue) : rule.range[0];
								rule.range[1]=(idvis.isNumber(rule.range[1]))? Math.max(rule.range[1],mappedValue) : rule.range[1];
							}
							else {
								rule.range=[mappedValue,mappedValue];
							}
						}
					}

					/* Point size range */
					if (CP.datasets[i].params.type==='point' && CP.datasets[i].params.pointSizeRule && !CP.datasets[i].params.pointSizeRule.userRange) { // => build range based on actual values 
						let rule=CP.datasets[i].params.pointSizeRule,
							mappedValue=newEl.propValues[rule.mapPropIndex];
						if (idvis.isNumber(mappedValue)) {
							if (rule.range) {
								rule.range[0]=Math.min(rule.range[0],mappedValue);
								rule.range[1]=Math.max(rule.range[1],mappedValue);
							}
							else {
								rule.range=[mappedValue,mappedValue];
							}
						}
					}

					if (datasetGroups[gIdx].layout==='stacked') baseValues[catIdx][gIdx]+=newEl.value;
					if (minCatPlotValue[axis]==undefined) { // assumes maxCatPlotValue is also undefined
						minCatPlotValue[axis]=newEl.value+newEl.baseValue;
						maxCatPlotValue[axis]=newEl.value+newEl.baseValue;
					}
					else {
						minCatPlotValue[axis]=Math.min(minCatPlotValue[axis],newEl.value+newEl.baseValue);
						maxCatPlotValue[axis]=Math.max(maxCatPlotValue[axis],newEl.value+newEl.baseValue);
					}
//console.log('min',minCatPlotValue[axis],'max',maxCatPlotValue[axis]);
				}
			}
		}
		// Adjusting pointSize scale(s) to display rounded numbers in legend
		if (this.scales && this.scales.pointSize && this.scales.pointSize.length) {
			this.scales.pointSize.forEach((rule) => {
				if (!rule.range || rule.transform) return; // declared wo range and never matched data
				let delta=rule.range[1]-rule.range[0];
				if (delta > 8 && delta % 8) { // use integer if possible
					const supl=8 - (delta % 8);
					delta+=supl;
					rule.range[1]+=supl;
				}
			});
		}
	}
	else { // histogram or area
		/* Find data range & compute bin width */
		let minValue=null,maxValue=null;
		for (let i=0; i<cpData.datasets.length; i++) {
			let set=cpData.datasets[i];
			for (let v=0; v<set.values.length; v++) {
//set.values[v]=-set.values[v];
				minValue=(minValue===null)? set.values[v] : Math.min(minValue,set.values[v]);
				maxValue=(maxValue===null)? set.values[v] : Math.max(maxValue,set.values[v]);
			}
		}
//console.log(minValue,maxValue);
        let numLinearbins=CP.categories.length;
        let correctMin=false, correctMax=false;
        if (CP.histLimits[0] !== null && CP.histLimits[0] > minValue) {
            minValue=CP.histLimits[0];
            numLinearbins--;
            correctMin=true;
        }
        if (CP.histLimits[1] !== null && CP.histLimits[1] < maxValue) {
            maxValue=CP.histLimits[1];
            numLinearbins--;
            correctMax=true;
        }
        if (correctMin || correctMax) {
            let tmpBinSize=(maxValue-minValue)/numLinearbins;
            if (correctMin) minValue-=tmpBinSize;
            if (correctMax) maxValue+=tmpBinSize;
        }
        
		/* Adjust (min/max)Value to optimized range */
		//histoClassOptRange=idvis.lib.getChartScaleRange(false,minValue,maxValue);
//console.log(histoClassOptRange);
		//minValue=histoClassOptRange[0];
		//maxValue=histoClassOptRange[2];
		
		let axis=(CP.horizontalCP)? 'Y' : 'X';
		minHistClassBinValue=minCatPlotValue[axis]=minValue;
		maxHistClassBinValue=maxCatPlotValue[axis]=maxValue;

		binSize=(maxValue-minValue)/CP.categories.length;
//console.log(minValue,maxValue,binSize);

		/* Distribute data in bins & build category elements */
		for (let i=0; i<cpData.datasets.length; i++) {
			let set=cpData.datasets[i];
			let axis=CP.datasets[i].params.axis;
			let gIdx=CP.datasets[i].params.group;
			for (let v=0; v<set.values.length; v++) {
				if (set.values[v]===undefined || set.values[v]===null) continue;
				//var catIdx=(correctMin && set.values[v] < CP.minThreshold)? 0 : (correctMax && set.values[v] > CP.maxThreshold)? CP.categories.length-1 : Math.floor((set.values[v]-minValue)/binSize); // =catIdx
				let val=(correctMin && set.values[v] < CP.histLimits[0])? minValue : (correctMax && set.values[v] > CP.histLimits[1])? maxValue : set.values[v];
                let catIdx=Math.floor((val-minValue)/binSize); // =catIdx
				catIdx=Math.min(catIdx,CP.categories.length-1); // highest value in list is max bin +1
				binCount[catIdx]++;
//if (set.values[v]<-12) {
//	console.log(set.values[v],catIdx);
//}
			}
			for (let catIdx=0; catIdx<binCount.length; catIdx++) {
				if (baseValues[catIdx][gIdx]===null) {baseValues[catIdx][gIdx]=0;}
				let newEl=new categoryElement(CP.datasets[i],catIdx,{value:binCount[catIdx],label:'Bin #'+(catIdx+1),id:catIdx},baseValues[catIdx][gIdx]);
				CP.datasets[i].data[catIdx]=newEl;
				if (datasetGroups[gIdx].layout==='stacked') baseValues[catIdx][gIdx]+=newEl.value;
				if (minCatPlotValue[axis]===undefined) { // assumes maxCatPlotValue is also undefined
					minCatPlotValue[axis]=newEl.value+newEl.baseValue;
					maxCatPlotValue[axis]=newEl.value+newEl.baseValue;
				}
				else {
					minCatPlotValue[axis]=Math.min(minCatPlotValue[axis],newEl.value+newEl.baseValue);
					maxCatPlotValue[axis]=Math.max(maxCatPlotValue[axis],newEl.value+newEl.baseValue);
				}
			}
		}
	}
	
	const axesUsedList=Object.keys(axesUsed);
    var toScale={};
	const noScale=function(value) {return value;};
	for (let i=0; i<axesUsedList.length; i++) {
		let axis=axesUsedList[i];
		if (axis.match('2')) { // Y2,X2
			//toScale[axis]=(cpData.toScale2)? cpData.toScale2 : function(value) {return value;};
			toScale[axis]=(CP.horizontalCP && (cpData.toScale || cpData.toScale1))? cpData.toScale || cpData.toScale1 : (CP.horizontalCP)? noScale : (cpData.toScale2)? cpData.toScale2 : noScale;
		}
		else { // Y,X
			//toScale[axis]=(cpData.toScale)? cpData.toScale : (cpData.toScale1)? cpData.toScale1 : function(value) {return value;};
			toScale[axis]=(CP.horizontalCP && cpData.toScale2)? cpData.toScale2 : (CP.horizontalCP)? noScale : (cpData.toScale || cpData.toScale1)? cpData.toScale || cpData.toScale1 : noScale;
		}
	}
//console.log(axesUsedList,toScale);
	this.outlierRule=(idvis.isNumber(cpData.outlierRule))? cpData.outlierRule : 'IQR';
	
	/********************* Data import *********************/
	/********** Add a feature ***********/
	this.addFeature=function(fData0) { // [axis,name,value,color,dashString,keepAbove1,roundValue]
		var fData={...fData0}; // shallow copy of fData0 (prevents modifying original object)
		var oldAxis=fData.axis || fData.axisX || 'N1';
		fData.axis=convertAxis(oldAxis);
		var feature=new idvis.feature(CP,fData); // [sh,minX*,maxX*,minY*,maxY*] *Only for custom & path
		if (feature.type==='line' && feature.editable) {
			this.editableThresholds.push(feature);
		}
//console.log(feature.name, feature.axis);
		//Add to list
		this.features.push(feature);
		//if (oldAxis.match('^N')) { // numerical axes
			if (minCatPlotValue[feature.axis]===undefined) {
				minCatPlotValue[feature.axis]=feature.startValue;
				maxCatPlotValue[feature.axis]=feature.startValue;
			}
			else {
				minCatPlotValue[feature.axis]=Math.min(minCatPlotValue[feature.axis],feature.startValue);
				maxCatPlotValue[feature.axis]=Math.max(maxCatPlotValue[feature.axis],feature.startValue);
			}
		//}
	};
	
	/********** Add a single box ***********/
	this.addBox=function(dsIdx,catIdx,boxInfo) {
		let gIdx=CP.datasets[dsIdx].params.group;
		if (baseValues[catIdx][gIdx]===null) {baseValues[catIdx][gIdx]=0;}
		let newBox=new box(CP.datasets[dsIdx],catIdx,boxInfo,baseValues[catIdx][gIdx]);
		CP.datasets[dsIdx].data[catIdx]=newBox;
		let axis=CP.datasets[dsIdx].params.axis;
		if (datasetGroups[gIdx].layout=='stacked') baseValues[catIdx][gIdx]+=newBox.median;
		if (minCatPlotValue[axis]===undefined) { // assumes maxCatPlotValue is also undefined
			minCatPlotValue[axis]=newBox.extremeValues[0]+newBox.baseValue;
			maxCatPlotValue[axis]=newBox.extremeValues[1]+newBox.baseValue;
		}
		else {
			minCatPlotValue[axis]=Math.min(minCatPlotValue[axis],newBox.extremeValues[0]+newBox.baseValue);
			maxCatPlotValue[axis]=Math.max(maxCatPlotValue[axis],newBox.extremeValues[1]+newBox.baseValue);
		}
	};
/*
	this.addData=function(dsIdx,data) { // OBSOLETE!!!!!!!!
		let gIdx=CP.datasets[dsIdx].params.group;
		if (CP.datasets[dsIdx].params.type=='histogram') {
			//code
		}
		else {
			for (var catIdx=0; catIdx<data.length; catIdx++) {
				if (!data[catIdx] || data[catIdx].length==0) continue;
				if (baseValues[catIdx][gIdx]==null) {baseValues[catIdx][gIdx]=0;}
				//var elData=data[catIdx].split(':');
				var newEl=new categoryElement(dsIdx,catIdx,data[catIdx],baseValues[catIdx][gIdx]);
				CP.datasets[dsIdx].data[catIdx]=newEl;
				var axis=CP.datasets[dsIdx].params.axis;
				if (datasetGroups[gIdx].layout=='stacked') baseValues[catIdx][gIdx]+=newEl.value;
				if (minCatPlotValue[axis]==undefined) { // assumes maxCatPlotValue is also undefined
					minCatPlotValue[axis]=newEl.value+newEl.baseValue;
					maxCatPlotValue[axis]=newEl.value+newEl.baseValue;
				}
				else {
					minCatPlotValue[axis]=Math.min(minCatPlotValue[axis],newEl.value+newEl.baseValue);
					maxCatPlotValue[axis]=Math.max(maxCatPlotValue[axis],newEl.value+newEl.baseValue);
				}
//console.log('min',minCatPlotValue[axis],'max',maxCatPlotValue[axis]);
			}
		}
	}
*/
    /*** Graphics variables ***/
	var topSpace=15; //(CP.datasets.length==1)? 20 : 30;
	var bottomSpace=20;
	var leftSpace=20;
	var rightSpace=20;
	var highlightZone=null;
	var minValue={},maxValue={},pix2valScale={}; //,thresholdValue,thresholdInputID;
	var canvasWidth,canvasHeight,catPlotX0,catPlotY0,catPlotW,catPlotH,catWidth,groupSize,catTitleFontSize,axisTitleFontSize,labelFontSize,rotateLabels=false;

	/********************* Box plot display (Raphael) *********************/
    this.draw=function() {
//console.log('Load='+(new Date().getTime()-t0));
//console.log(minCatPlotValue,maxCatPlotValue);
		/* DIVs */
		/*
		if (this.editThreshold || this.datasets.length > 1) {
			canvasDivID=this.divID+'_canvas';
			formDivID=this.divID+'_form';
			this.formDiv=document.getElementById(formDivID);
		}
		else {
			canvasDivID=this.divID;
		}
		*/

		/* Computing histoClass category axis scale (Must be done before numerical axes) */
		let rangeBeforeBins=0,rangeAfterBins=0;
		//let [axis,categorySize]=(CP.horizontalCP)? ['Y',catPlotH] : ['X',catPlotW];
		let axis=(CP.horizontalCP)? 'Y' : 'X';
		if (histoClass) {
			let optimize=(minCatPlotValue[axis] < minHistClassBinValue || maxCatPlotValue[axis] > maxHistClassBinValue)? true : false;
			histoClassOptRange=idvis.lib.getChartScaleRange(optimize,minCatPlotValue[axis],maxCatPlotValue[axis]);
			minValue[axis]=histoClassOptRange[0];
			maxValue[axis]=histoClassOptRange[2];
			rangeBeforeBins=Math.max(0,minHistClassBinValue-minValue[axis]);
			rangeAfterBins=Math.max(0,maxValue[axis]-maxHistClassBinValue);
		}
		/* Computing Graphic area size */
		const paddingFactor=(cpData.padding===undefined)? 1.2 : (cpData.padding >= 1)? 100/(100-cpData.padding) : 1/(1-cpData.padding); // >=1 assume % otherwise fraction
		const catElemSize=cpData.boxSize || cpData.barSize || cpData.binSize; // can be undefined
//console.log('paddingFactor',paddingFactor);
		if (CP.horizontalCP) {
			if (catElemSize) {
				groupSize=Math.max(3,catElemSize);
				catWidth=Math.round(groupSize*paddingFactor*datasetGroups.length);
				spaceBeforeBins=(rangeBeforeBins)? (groupSize*paddingFactor/binSize)*rangeBeforeBins : 0;
				spaceAfterBins=(rangeAfterBins)? (groupSize*paddingFactor/binSize)*rangeAfterBins : 0;
				catPlotH=(catWidth*this.categories.length) + spaceBeforeBins + spaceAfterBins;
				pix2valScale[axis]=(maxValue[axis]-minValue[axis])/catPlotH;
			}
			else if (cpData.height) {
				catPlotH=cpData.height;
				pix2valScale[axis]=(maxValue[axis]-minValue[axis])/catPlotH;
				spaceBeforeBins=(rangeBeforeBins)? rangeBeforeBins/pix2valScale[axis] : 0;
				spaceAfterBins=(rangeAfterBins)? rangeAfterBins/pix2valScale[axis] : 0;	
				catWidth=(catPlotH-spaceBeforeBins-spaceAfterBins)/this.categories.length;
				groupSize=Math.round(catWidth/(paddingFactor*datasetGroups.length));
			}
			else { // auto
				catWidth=Math.min(150,Math.max(11,750/this.categories.length));
				groupSize=Math.round(catWidth/(paddingFactor*datasetGroups.length));
				spaceBeforeBins=(rangeBeforeBins)? (groupSize*paddingFactor/binSize)*rangeBeforeBins : 0;
				spaceAfterBins=(rangeAfterBins)? (groupSize*paddingFactor/binSize)*rangeAfterBins : 0;
				catPlotH=(catWidth*this.categories.length) + spaceBeforeBins + spaceAfterBins;
				pix2valScale[axis]=(maxValue[axis]-minValue[axis])/catPlotH;	
			}
			catPlotW=(cpData.width)? cpData.width : 250;
			axisTitleFontSize=(catPlotW > 250)? 14 : 10;
			labelFontSize=(catWidth > 20)? 14 : 10;
			let labelSpaceX;
			if (histoClass) {
				labelSpaceX=40;
			}
			else {
				let tmpCanvas=Raphael(mainDivID,0,0),
					ml=tmpCanvas.text(0,0,maxLabel).attr({'font-size':labelFontSize,'font-weight':'bold'}),
					maxLabelWidth=Math.round(ml.getBBox().width);
				tmpCanvas.remove();
				labelSpaceX=maxLabelWidth + 10;
			}
			catPlotX0=leftSpace+labelSpaceX+1;
			catPlotY0=topSpace+0.5;
			if (axesUsed.X2) catPlotY0+=40;
			canvasHeight=catPlotY0+catPlotH+bottomSpace;
			if (axesUsed.X) canvasHeight+=30;
			canvasWidth=leftSpace+labelSpaceX+catPlotW+rightSpace;
		}
		else { // vertical (normal case)
			if (catElemSize) {
				groupSize=Math.max(3,catElemSize);
				catWidth=groupSize*paddingFactor*datasetGroups.length;
				spaceBeforeBins=(rangeBeforeBins)? (groupSize*paddingFactor/binSize)*rangeBeforeBins : 0;
				spaceAfterBins=(rangeAfterBins)? (groupSize*paddingFactor/binSize)*rangeAfterBins : 0;
				catPlotW=(catWidth*this.categories.length) + spaceBeforeBins + spaceAfterBins;
				pix2valScale[axis]=(maxValue[axis]-minValue[axis])/catPlotW;
			}
			else if (cpData.width) {
				catPlotW=cpData.width;
				pix2valScale[axis]=(maxValue[axis]-minValue[axis])/catPlotW;
				spaceBeforeBins=(rangeBeforeBins)? rangeBeforeBins/pix2valScale[axis] : 0;
				spaceAfterBins=(rangeAfterBins)? rangeAfterBins/pix2valScale[axis] : 0;
				catWidth=(catPlotW-spaceBeforeBins-spaceAfterBins)/this.categories.length;
				groupSize=(catWidth/(paddingFactor*datasetGroups.length));
			}
			else { // auto
				catWidth=Math.min(150,Math.max(11,750/this.categories.length));
				groupSize=catWidth/(paddingFactor*datasetGroups.length);
				spaceBeforeBins=(rangeBeforeBins)? (groupSize*paddingFactor/binSize)*rangeBeforeBins : 0;
				spaceAfterBins=(rangeAfterBins)? (groupSize*paddingFactor/binSize)*rangeAfterBins : 0;
				catPlotW=(catWidth*this.categories.length) + spaceBeforeBins + spaceAfterBins;
				pix2valScale[axis]=(maxValue[axis]-minValue[axis])/catPlotW;
			}
			catPlotH=(cpData.height)? cpData.height : 250;
			axisTitleFontSize=(catPlotW > 250)? 14 : 10;
			labelFontSize=(catWidth > 20)? 14 : 10;
			let labelSpaceY; // space for category title and label!!!
			if (histoClass) {
				labelSpaceY=25;
			}
			else {
				let tmpCanvas=Raphael(mainDivID,0,0),
				    ml=tmpCanvas.text(0,0,maxLabel).attr({'font-size':labelFontSize,'font-weight':'bold'}),
				    maxLabelWidth=Math.round(ml.getBBox().width);
				tmpCanvas.remove();
				if (maxLabelWidth > catWidth-5) {
					rotateLabels=true;
					labelSpaceY=maxLabelWidth+10;
				}
				else {
					labelSpaceY=35;
				}
			}
			catPlotX0=leftSpace+0.5;
			if (axesUsed.Y) catPlotX0+=50; //+Math.max(30,labelSpaceY)+1;
			catPlotY0=topSpace+1.5;
			canvasWidth=catPlotX0-1+catPlotW+rightSpace;
			if (axesUsed.Y2) canvasWidth+=50;
			canvasHeight=topSpace+catPlotH+labelSpaceY+bottomSpace;
		}

		catTitleFontSize=(catWidth*this.categories.length > 250)? 18 : 14;
		const maxPointDiameter=Math.min(0.75 * groupSize, 50);

		/* Compute space required for scale legends */
		let legendStartX=canvasWidth;
		if (this.showScales && typeof(this.scales)==='object') {
			let maxLabel='';
			['color','pointSize'].forEach((prop) => {
				if (this.scales[prop] && this.scales[prop].length) {
					this.scales[prop].forEach((rule) => {
						if (!rule.range) return; // declared wo range and never matched data
						let ruleLabel=(rule.transform)? rule.transform.label : rule.map;
						if (ruleLabel.length > maxLabel.length) maxLabel=ruleLabel;
					});
				}
			});
			let tmpCanvas=Raphael(mainDivID,0,0),
				ml=tmpCanvas.text(0,0,maxLabel).attr({'font-size':14,'font-weight':'bold'});
			canvasWidth+=Math.max(ml.getBBox().width,groupSize)+20; // groupSize <=> biggest point to plot
			tmpCanvas.remove();
		}

		/* Needed by idvis in case editable thresholds */
		CP.chartSettings.current[axis]={
			minRange: minValue[axis],
			//maxRange: maxValue[axis],
			pix2valRatio: pix2valScale[axis]
		};
		/* DIVs for form and plot*/
		let formPosition=(cpData.formPosition && cpData.formPosition.toLowerCase().match('^(t|b|r|l)'))? cpData.formPosition.toLowerCase() : (CP.editableThresholds.length || this.hasBoxes || CP.hasPoints > 1)? 'right' : 'top';
		idvis.lib.drawChartLayout(CP,formPosition,canvasWidth,canvasDivID,formDivID);
		
		/*** Form menu ***/
		CP.formDiv=document.getElementById(formDivID);
		var customFormElements=[];
//		customFormElements[10]=[false]; // hide idvis default Dataset selection menu
		/* Custom Dataset selection */
		/*
		if (CP.hasPoints > 1) {
			let htmlString='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Points:</B></LEGEND>';
			let dsPointCount=0;
			for (let i=0; i<CP.datasets.length; i++) {
				if (CP.datasets[i].params.type != 'point') continue;
				if (dsPointCount > 0) {htmlString+='<BR>';}
				htmlString+='<INPUT type="checkbox" value="'+i+'" onclick="idvis.lib.setDatasetDisplay(idvis.registeredCharts['+CP.chartID+'],this)" checked><FONT style="color:'+CP.datasets[i].params.color+'">'+CP.datasets[i].params.name+'</FONT>';
				dsPointCount++;
			}
			htmlString+="</FIELDSET>\n";
			customFormElements[1]=htmlString; //[0,10,20....] reserved to control default idvis form elements
		}
		*/
		/* Show/hide box points */
		if (CP.hasBoxes) {
			htmlString='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Boxes:</B></LEGEND>';
			htmlString+='<INPUT type="checkbox" value="1" onclick="idvis.registeredCharts['+CP.chartID+'].showBoxValues(this.checked)"><B>Show all values</B>';
			htmlString+="</FIELDSET>\n";
			customFormElements[21]=htmlString;		
		}
		idvis.lib.initializeForm(CP,customFormElements); // form elements creation

		/* Canvas chart */
		this.canvas=Raphael(canvasDivID,canvasWidth,canvasHeight);

		/* Back panel */
		this.canvas.rect(0,0,canvasWidth,canvasHeight,0).attr({fill:'#fff',stroke:'#000'});

		/* Chart area */
		this.plotArea=this.canvas.rect(catPlotX0,catPlotY0,catPlotW,catPlotH,0).attr({fill:'#F3F3F3',stroke:'#000'}); //,cursor:"crosshair"

		/* highlight zone */
		highlightZone=(CP.horizontalCP)? this.canvas.rect(catPlotX0,catPlotY0,catPlotW,catWidth,0) : this.canvas.rect(catPlotX0,catPlotY0,catWidth,catPlotH,0);
		highlightZone.attr({fill:'#F00',stroke:'none','fill-opacity':0.1}).hide();

		/* Dataset legends (skip for point if > 1) */
/* OBSOLETE
		if (CP.datasets.length > 1) {
			let ly=17, lx=catPlotX0;
			for (let i=0; i < CP.datasets.length; i++) {
				if (CP.hasPoints > 1 && CP.datasets[i].params.type == 'point') continue;
				// Dataset type logo
				drawIcon(CP.datasets[i].params.type,lx,ly-8.5,CP.datasets[i].params.color);
				// Dataset name
				lx+=20;
				let lt=CP.canvas.text(lx,ly,CP.datasets[i].params.name).attr({'font-size':labelFontSize,'font-weight':'bold','text-anchor':'start',fill:CP.datasets[i].params.color});
				lx+=lt.getBBox().width+20;
			}
		}
*/
		/* Check toScale function(s) */
		for (let i=0; i<axesUsedList.length; i++) {
			let axis=axesUsedList[i];
			let failed=false,failedValue=null;
			if (!idvis.isNumber(toScale[axis](minCatPlotValue[axis]))) {failed=true; failedValue=minCatPlotValue[axis];}
			else if (!idvis.isNumber(toScale[axis](maxCatPlotValue[axis]))) {failed=true; failedValue=maxCatPlotValue[axis];}
			if (failed) {
				//var axisNumber=(axis.match('2'))? '2' : '1';
				alert('ERROR: Conversion of "'+failedValue+'" for axis "'+axis+'" is not a valid number.');
				this.canvas.text(catPlotX0+(catPlotW/2),catPlotY0+(catPlotH/2),'Drawing has failed.').attr({'font-weight':'bold','font-size':14,fill:'#DD0000'});
				return;
			}
		}
		
		/* Compute axes scale */
		var numericalSize=(CP.horizontalCP)? catPlotW : catPlotH;
		var optMinV={},tickSizeV={};
		var hasNegativeValues=false;
		for (let i=0; i<axesUsedList.length; i++) { // Numerical axi(e)s
			let axis=axesUsedList[i],
				toScaleMinVal=(startAt0)? Math.min(0,toScale[axis](minCatPlotValue[axis])) : toScale[axis](minCatPlotValue[axis]),
				optRange=idvis.lib.getChartScaleRange(true,toScaleMinVal,toScale[axis](maxCatPlotValue[axis]),toScale[axis](numericalSize)); // requires chartLibrary
			minValue[axis]=(startAt0 && toScaleMinVal >= 0)? 0 : optRange[0];
			if (toScaleMinVal < 0) hasNegativeValues=true;
			optMinV[axis]=optRange[1];
			maxValue[axis]=optRange[2];
			tickSizeV[axis]=optRange[3];
			pix2valScale[axis]=(maxValue[axis]-minValue[axis])/toScale[axis](numericalSize);
			numberOfAxes++;
			/* Needed by idvis in case editable thresholds */
			CP.chartSettings.current[axis]={
				minRange: minValue[axis],
				//maxRange: maxValue[axis],
				pix2valRatio: pix2valScale[axis]
			};
		}
		
		/* Adjust 1 axis scale so both 0 lines overlap TODO: adjust axies ticks also */
		if (numberOfAxes==2 && hasNegativeValues) { // Negative values
			let minValProp={};
			for (let i=0; i<axesUsedList.length; i++) {
				let axis=axesUsedList[i];
				minValProp[axis]=(0 - minValue[axis]) / (maxValue[axis] - minValue[axis]);
			}
			let axis1,axis2;
			if (CP.horizontalCP) {axis1='X'; axis2='X2';}
			else {axis1='Y'; axis2='Y2';}
			if (minValProp[axis1] > minValProp[axis2]) { // adjust axis2
				minValue[axis2]=minValue[axis1]*maxValue[axis2]/maxValue[axis1];
				pix2valScale[axis2]=(maxValue[axis2]-minValue[axis2])/toScale[axis2](numericalSize);
				CP.chartSettings.current[axis2]={ // overwrite previous
					minRange: minValue[axis2],
					pix2valRatio: pix2valScale[axis2]
				};
//console.log(axis2,(0 - minValue[axis2]) / (maxValue[axis2] - minValue[axis2]))
			}
			else { // adjust axis1
				minValue[axis1]=minValue[axis2]*maxValue[axis1]/maxValue[axis2];
				pix2valScale[axis1]=(maxValue[axis1]-minValue[axis1])/toScale[axis1](numericalSize);
				CP.chartSettings.current[axis1]={ // overwrite previous
					minRange: minValue[axis1],
					pix2valRatio: pix2valScale[axis1]
				};
//console.log(axis1,(0 - minValue[axis1]) / (maxValue[axis1] - minValue[axis1]))
			}
		}

		var catPlotStartY0=catPlotY0+catPlotH;
		var zeroValuePos;

		/** Axes titles & ticks **/
		/* Titles */
		//var axisXtext,axisYtext;
		var axisText={}, axisFontSize={};
		if (CP.horizontalCP) {
			//axisXtext=(cpData.valueAxes && cpData.valueAxes.title)? cpData.valueAxes.title : 'Values';
			//axisYtext=(cpData.categories.title)? cpData.categories.title : 'Categories';
			for (let i=0; i<axesUsedList.length; i++) { //X,X2
				let axis=axesUsedList[i],
					shortTitle;
				if (cpData.valueAxes) {
					if (axis==='X') {
						axisText[axis]=cpData.valueAxes.title || cpData.valueAxes.title1;
						shortTitle=cpData.valueAxes.shortTitle || cpData.valueAxes.shortTitle1; // can undef
					}
					else { // X2
						axisText[axis]=cpData.valueAxes.title2;
						shortTitle=cpData.valueAxes.shortTitle2;
					}
				}
				if (!axisText[axis]) axisText[axis]='Values';
				CP.axesTitles[axis]={title:axisText[axis]};
				if (shortTitle) CP.axesTitles[axis].shortTitle=shortTitle;
			}
			axisText.Y=cpData.categories.title || 'Categories';
			axisFontSize.Y=catTitleFontSize;
		}
		else { // default
			for (let i=0; i<axesUsedList.length; i++) {  // Y,Y2
				let axis=axesUsedList[i],
					shortTitle;
				if (cpData.valueAxes) {
					if (axis==='Y') {
						axisText[axis]=cpData.valueAxes.title || cpData.valueAxes.title1;
						shortTitle=cpData.valueAxes.shortTitle || cpData.valueAxes.shortTitle1; // can undef
					}
					else { // Y2
						axisText[axis]=cpData.valueAxes.title2;
						shortTitle=cpData.valueAxes.shortTitle2;
					}
				}
				if (!axisText[axis]) axisText[axis]='Values';
				CP.axesTitles[axis]={title:axisText[axis]};
				if (shortTitle) CP.axesTitles[axis].shortTitle=shortTitle;
			}
			axisText.X=cpData.categories.title || 'Categories';
			axisFontSize.X=catTitleFontSize;
		}
		//if (axisText.X) this.canvas.text(catPlotX0+catPlotW/2,catPlotY0-30,axisText.X).attr({'font-weight':'bold','font-size':14});
		//let posYTextX2=(CP.horizontalCP)? catPlotY0+catPlotH+30 : canvasHeight-15;
		//if (axisText.X2) this.canvas.text(catPlotX0+catPlotW/2,posYTextX2,axisText.X2).attr({'font-weight':'bold','font-size':14});
		if (axisText.X) this.canvas.text(catPlotX0+catPlotW/2,canvasHeight-15,axisText.X).attr({'font-weight':'bold','font-size':axisFontSize.X || axisTitleFontSize});		
		if (axisText.X2) this.canvas.text(catPlotX0+catPlotW/2,catPlotY0-30,axisText.X2).attr({'font-weight':'bold','font-size':axisTitleFontSize});		
		let ty=catPlotY0+(catPlotH/2);
		let rotateY=(CP.horizontalCP)? 90 : -90;
		if (axisText.Y) this.canvas.text(15,ty,axisText.Y).attr({'font-weight':'bold','font-size':axisFontSize.Y || axisTitleFontSize}).rotate(rotateY,15,ty);
		if (axisText.Y2) this.canvas.text(canvasWidth-15,ty,axisText.Y2).attr({'font-weight':'bold','font-size':axisTitleFontSize}).rotate(-90,canvasWidth-15,ty);
		
		/* ticks */
		var histoTicks=false;
		for (let i=0; i<axesUsedList.length; i++) {
			let axis=axesUsedList[i],
				posX,posY;
			//var fixed=(tickSizeV[axis] < 1)? (Math.round(1/tickSizeV[axis])+'').length : (Math.round(tickSizeV[axis])==tickSizeV[axis])? 0 : 1;
			if (CP.horizontalCP) { // X ticks
				let textY;
				if (axis==='X2') {
					posY=catPlotY0-6;
					textY=posY-16;
				}
				else {
					posY=catPlotY0+catPlotH;
					textY=posY+1;
				}
				let tick=optMinV[axis];
				while (tick <= maxValue[axis]) {
					posX=catPlotX0+Math.round((tick-minValue[axis])/pix2valScale[axis]);
					if (posX >= catPlotX0) {
						this.canvas.path('M'+posX+' '+posY+' l0 5');
						this.canvas.text(posX,textY+10,idvis.formatNumber(tick));
					}
					tick+=tickSizeV[axis];
				}
				/* 0 line for boxplot */
				if (startAt0 && hasNegativeValues) {
					if (!zeroValuePos) { // draw only once
						zeroValuePos=catPlotX0-Math.round(minValue[axis]/pix2valScale[axis]);
						this.canvas.path('M'+zeroValuePos+' '+(catPlotStartY0-catPlotH)+'l0 '+catPlotH);
					}
				}
				else{zeroValuePos=catPlotX0;}
			}
			else { // Y axis (default)
				let textX,textAnch;
				if (axis==='Y') {
					posX=catPlotX0-6;
					textX=posX-1;
					textAnch='end';
				}
				else {
					posX=catPlotX0+catPlotW;
					textX=posX+7;
					textAnch='start';
				}
				let tick=optMinV[axis];
				while (tick <= maxValue[axis]) {
					posY=catPlotStartY0-Math.round((tick-minValue[axis])/pix2valScale[axis]);
					if (posY <= catPlotStartY0) {
						this.canvas.path('M'+posX+' '+posY+' l5 0');
						this.canvas.text(textX,posY,idvis.formatNumber(tick)).attr({'text-anchor':textAnch});
					}
					tick+=tickSizeV[axis];
				}
				/* 0 line for boxplot */
				if (startAt0 && hasNegativeValues) {
					if (!zeroValuePos) { // draw only once
						zeroValuePos=catPlotStartY0+Math.round(minValue[axis]/pix2valScale[axis]);
						this.canvas.path('M'+catPlotX0+' '+zeroValuePos+'l'+catPlotW+' 0');
					}
				}
				else{zeroValuePos=catPlotStartY0;}
			}
			if (histoClass && !histoTicks) { // do only once
				histoTicks=true;
				let [minValueHC,optMinHC,maxValueHC,tickSizeHC]=histoClassOptRange;
				//var posX,posY;
				//fixed=(tickSizeHC < 1)? (Math.round(1/tickSizeHC)+'').length : (Math.round(tickSizeHC)==tickSizeHC)? 0 : 1;
				if (CP.horizontalCP) { // Y axis ticks
					let pix2valScaleHC=(maxValueHC-minValueHC)/catPlotH;
					posX=catPlotX0-6;
					let tick=optMinHC;
					while (tick <= maxValueHC) {
						posY=catPlotY0+Math.round((tick-minValueHC)/pix2valScaleHC);
						if (posY >= catPlotY0) {
							this.canvas.path('M'+posX+' '+posY+' l5 0');
							//this.canvas.text(posX-1,posY,tick.toFixed(fixed)).attr({'text-anchor':'end'});
							this.canvas.text(posX-1,posY,idvis.formatNumber(tick)).attr({'text-anchor':'end'});
							
						}
						tick+=tickSizeHC;
					}
				}
				else { // X axis ticks
					let pix2valScaleHC=(maxValueHC-minValueHC)/catPlotW;
					posY=catPlotStartY0+1;
					posX=catPlotX0;
					let tick=optMinHC;
					while (tick <= maxValueHC) {
						posX=catPlotX0+Math.round((tick-minValueHC)/pix2valScaleHC);
						if (posX >= catPlotX0) {
							this.canvas.path('M'+posX+' '+posY+' l0 5');
							//this.canvas.text(posX,posY+10,tick.toFixed(fixed));
							this.canvas.text(posX,posY+10,idvis.formatNumber(tick));
						}
						tick+=tickSizeHC;
					}
				}	
			}
		}

		/***** Drawing boxes/points & labels *****/
		var catInnerWidth=groupSize*datasetGroups.length,
		    boxSize=(datasetGroups.length > 1)? groupSize/paddingFactor : groupSize,
		    halfBoxSize=boxSize/2,
		    whiskersW=Math.max(7,Math.round(halfBoxSize)), // boxplot only
		    catCenterPos=(CP.horizontalCP)? catPlotY0+spaceBeforeBins-(catWidth/2) : catPlotX0+spaceBeforeBins-(catWidth/2),
            zig = new idvis.ziggurat(); // for normally distributed random numbers
		for (let i=0; i<CP.categories.length; i++) {
			catCenterPos+=catWidth;
			
			/* category ticks */
			if (!histoClass) {
				//if (i > 0) {
				if (CP.horizontalCP) {
					let tX=catPlotX0-5,
					    tY=catPlotY0+(catWidth*i);
					CP.canvas.path('M'+tX+' '+tY+' l10 0');
					if (i==CP.categories.length-1) { // add last tick
						CP.canvas.path('M'+tX+' '+(catPlotY0+catPlotH)+' l10 0');
					}
				}
				else {
					let tX=catPlotX0+(catWidth*i),
					    tY=catPlotY0+catPlotH;
					CP.canvas.path('M'+tX+' '+tY+' l0 5');
					if (i==CP.categories.length-1) { // add last tick
						CP.canvas.path('M'+(catPlotX0+catPlotW)+' '+tY+' l0 5');
					}
				}
				//}

				/* labels */
				let lx,ly,la;
				if (CP.horizontalCP) {
					lx=catPlotX0-5;
					ly=catCenterPos;
					la=0;
				}
				else {
					lx=catCenterPos;
					ly=(rotateLabels)? catPlotStartY0+5 : catPlotStartY0+15;
					la=-90; //-45;
				}
				CP.categories[i].labelSvg=CP.canvas.text(lx,ly,CP.categories[i].label)
				.attr({'font-size':labelFontSize,'font-weight':'bold','fill':CP.categories[i].labelColor})
				.data('js',CP.categories[i])
				.hover(function(e){setCatEmphasis(this,'on',e)},function(){setCatEmphasis(this,'off')})
				.click(function(){categoryOnClick(this.data('js').externalID)});
				if (rotateLabels || CP.horizontalCP) {CP.categories[i].labelSvg.attr({'text-anchor':'end'})}
				if (rotateLabels) {CP.categories[i].labelSvg.attr({'text-anchor':'end'}).rotate(la,lx,ly)}
				//if (boxOnClick) {CP.categories[i].labelSvg.click(function() {boxOnClick(this.data('js').externalID)});}
			}
			
			var groupCenterPos=catCenterPos-(catInnerWidth/2)-(groupSize/2);
			for (let g=0; g<datasetGroups.length; g++) {
				groupCenterPos+=groupSize;
				for (let d=0; d<datasetGroups[g].content.length; d++) {
					let dsIdx=datasetGroups[g].content[d],
					    axis=datasetGroups[g].axis || CP.datasets[dsIdx].params.axis,  // only defined for 'stacked'
					    catElem=CP.datasets[dsIdx].data[i];
					if (!catElem) continue; // no dataset element in this category
					let baseShift=Math.round(catElem.baseValue/pix2valScale[axis]); // 0 if no stacking;
					
					/* Color calculation */					
					if (CP.datasets[dsIdx].params.colorRule) {
						let rule=CP.datasets[dsIdx].params.colorRule,
							transformFunc=(rule.transform)? rule.transform.func : (x)=>x, // arrow function
							transR0=transformFunc(rule.range[0]),
							transR1=transformFunc(rule.range[1]),
							transRange=(transR0 <= transR1)? [transR0,transR1] : [transR1,transR0],
							delta=transRange[1]-transRange[0],
							mappedValue=(catElem.propValues && idvis.isNumber(catElem.propValues[rule.mapPropIndex]))? catElem.propValues[rule.mapPropIndex] : rule.range[0];
						mappedValue=Math.min(transRange[1],Math.max(transformFunc(mappedValue),transRange[0])); // make sure value is within range

						// let relValueDown=Math.round(70*(Math.max(0,Math.min(delta,rule.range[1]-mappedValue)))/delta),
						// 	relValueUp=Math.round(100*(Math.max(0,Math.min(delta,mappedValue-rule.range[0])))/delta);
						// catElem.color = 'rgb('+relValueDown+'%,'+relValueDown+'%,'+relValueUp+'%)'; // 70%,70%,0% <-> 0%,0%,100% (dark yellow <-> blue)
						let relValueDown=Math.max(0,Math.min(delta,transRange[1]-mappedValue))/delta,
							relValueUp=Math.max(0,Math.min(delta,mappedValue-transRange[0]))/delta;
						catElem.color = scaleGradients[rule.gradient](relValueDown,relValueUp);
					} // else catElem.color already defined at catElem initialization

					/** box **/
					if (CP.datasets[dsIdx].params.type==='box') {
						let box,median,whiskers,
						    lowerWhSize=Math.round((toScale[axis](catElem.lowerQuartile)-toScale[axis](catElem.lowerWhisker))/pix2valScale[axis]),
						    upperWhSize=Math.round((toScale[axis](catElem.upperWhisker)-toScale[axis](catElem.upperQuartile))/pix2valScale[axis]);
						if (CP.horizontalCP) {
							let bx=catPlotX0+Math.round((toScale[axis](catElem.lowerQuartile)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])+baseShift,
							    by=Math.round(groupCenterPos-halfBoxSize),
							    bw=Math.round((toScale[axis](catElem.upperQuartile)-toScale[axis](catElem.lowerQuartile))/pix2valScale[axis])+baseShift;
							box=CP.canvas.rect(bx,by,bw,boxSize);
							let mx=catPlotX0+Math.round((toScale[axis](catElem.median)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])+baseShift;
							median=CP.canvas.path('M'+mx+' '+by+' l0 '+boxSize);
							let wlx=catPlotX0+Math.round((toScale[axis](catElem.lowerWhisker)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])+baseShift,
							    wy=Math.round(groupCenterPos-(whiskersW/2)),
							    wux=catPlotX0+Math.round((toScale[axis](catElem.upperWhisker)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])+baseShift;
							whiskers=CP.canvas.path('M'+wlx+' '+wy+' l0 '+whiskersW+' M'+wlx+' '+groupCenterPos+' l'+lowerWhSize+' 0 M'+wux+' '+wy+' l 0 '+whiskersW+' M'+wux+' '+groupCenterPos+' l-'+upperWhSize+' 0');
							catElem.x=wux;
							catElem.y=groupCenterPos;
							let meax=catPlotX0+Math.round((toScale[axis](catElem.mean)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])+baseShift;
							catElem.meanSvg=CP.canvas.path('M'+meax+' '+by+' l0 '+boxSize);
						}
						else { // default
							var bx=Math.round(groupCenterPos-halfBoxSize),
							    by=catPlotStartY0-Math.round((toScale[axis](catElem.upperQuartile)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])-baseShift,
							    bh=Math.round((toScale[axis](catElem.upperQuartile)-toScale[axis](catElem.lowerQuartile))/pix2valScale[axis]);
							box=this.canvas.rect(bx,by,boxSize,bh);
							let my=catPlotStartY0-Math.round((toScale[axis](catElem.median)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])-baseShift;
							median=CP.canvas.path('M'+bx+' '+my+' l'+boxSize+' 0');
							let wx=Math.round(groupCenterPos-(whiskersW/2)),
							    wly=catPlotStartY0-Math.round((toScale[axis](catElem.lowerWhisker)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])-baseShift,
							    wuy=catPlotStartY0-Math.round((toScale[axis](catElem.upperWhisker)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])-baseShift;
							whiskers=CP.canvas.path('M'+wx+' '+wly+' l'+whiskersW+' 0 M'+groupCenterPos+' '+wly+' l0 -'+lowerWhSize+' M'+wx+' '+wuy+' l'+whiskersW+' 0 M'+groupCenterPos+' '+wuy+' l0 '+upperWhSize);					
							catElem.x=groupCenterPos;
							catElem.y=wuy;
							let meay=catPlotStartY0-Math.round((toScale[axis](catElem.mean)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])-baseShift;
							catElem.meanSvg=CP.canvas.path('M'+bx+' '+meay+' l'+boxSize+' 0');
						}
						//box.attr({fill:catElem.color,stroke:catElem.color,'fill-opacity':0.3});
						median.attr({'stroke-width':3}); //stroke:catElem.color,
						catElem.meanSvg.attr({'stroke-width':3,'stroke-dasharray':'.',stroke:catElem.color});
						//whiskers.attr({stroke:catElem.color,'stroke-width':2});
						catElem.svg=CP.canvas.set(box,median,whiskers).attr({fill:catElem.color,stroke:catElem.color,'fill-opacity':CP.datasets[dsIdx].params.opacity})
						.data('js',catElem)
						.hover(function(e){setBoxEmphasis(this,'on',e)},function(){setBoxEmphasis(this,'off')}) // this=any element in set (not set!)
						.click(function() {boxElemClicked(this.data('js'))});
						if (catElem.excluded) {catElem.svg.attr({'stroke-dasharray':'- ','stroke-opacity':0.5});}
						//if (boxOnClick) {catElem.svg.click(function() {boxOnClick(this.data('js').externalID)});}
			
						/* box points with outliers */
                        let minPos=groupCenterPos-halfBoxSize,
                            maxPos=minPos+boxSize,
                            rangeUp=catElem.upperWhisker-catElem.upperQuartile,
                            rangeLow=catElem.lowerQuartile-catElem.lowerWhisker;
						for (let j=0; j< catElem.pointList.length; j++) {
							let point=catElem.pointList[j],
                                range=(point.isOutlier)? 0 : (point.isExcluded)? 0.25 : (point.value > catElem.upperQuartile)? (catElem.upperWhisker-point.value)/rangeUp : (point.value < catElem.lowerQuartile)? (point.value-catElem.lowerWhisker)/rangeLow : 1,
                                pos=Math.round( Math.max( minPos,Math.min( maxPos,groupCenterPos + zig.nextGaussian()*halfBoxSize*range/3 ) ) ); // 3 -> +/-3 std dev = 99.7% of data
							if (CP.horizontalCP) {
								let px=catPlotX0+Math.round((toScale[axis](point.value)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])+baseShift;
								point.svg=point.point=CP.canvas.circle(px,pos,3);
							}
							else {
								let py=catPlotStartY0-Math.round((toScale[axis](point.value)-minValue[axis]+catElem.baseValue)/pix2valScale[axis])-baseShift;
                                point.svg=point.point=CP.canvas.circle(pos,py,3);
							}
							let pColor=(point.color)? point.color : catElem.color; // boxplot excluded points have exclusion color
							point.point.attr({fill:pColor,stroke:pColor,'fill-opacity':0.4})
							.data({js:point,showLabel:null})
							//.hover(function(){setBoxPointEmphasis(this.data('js'),'on')},function(){setBoxPointEmphasis(this.data('js'),'off')})
							//.click(function(){selectBoxPoint(this.data('js'))}); //,'auto'
							.hover(function(e){let type=(idvis.modKeyPressed)? 'max' : 'min'; idvis.lib.setLabelDisplay(this,'on',type,e);},function(){idvis.lib.setLabelDisplay(this,'off');})
							.click(function(){idvis.lib.setPointSelection(this,'auto','max');});
							if (catElem.excluded) {point.svg.attr({'stroke-dasharray':'. ','stroke-opacity':0.5});}
                            if (point.isOutlier===false && point.isExcluded===false) {point.svg.hide();}
						}
					}
					
					/** Bar or Histogram **/
					else if (CP.datasets[dsIdx].params.type.match('bar|histo')) {
						//var bh=Math.round((toScale[axis](catElem.value)-minValue[axis])/pix2valScale[axis]);
						let bValue=toScale[axis](catElem.value),
						    bh=Math.round(bValue/pix2valScale[axis]),
						    bar;
						if (CP.horizontalCP) {
							let by=Math.round(groupCenterPos-halfBoxSize);
							if (bValue >= 0) {
								//var bx=catPlotX0+baseShift;
								let bx=zeroValuePos+baseShift;
								bar=CP.canvas.rect(bx+0.5,by,bh-0.5,boxSize); // +/-0.5 not to hide chart border
							}
							else {
								let bx=zeroValuePos+bh+baseShift;
								bar=CP.canvas.rect(bx+0.5,by,0.5-bh,boxSize); // +/-0.5 not to hide chart border
							}
						}
						else { // default
							let bx=Math.round(groupCenterPos-halfBoxSize);
							if (bValue >= 0) {
								let by=zeroValuePos-bh-baseShift;
								bar=CP.canvas.rect(bx,by-0.5,boxSize,bh-0.5); // +/-0.5 not to hide chart border
							}
							else {
								let by=zeroValuePos-baseShift;
								bar=CP.canvas.rect(bx,by-0.5,boxSize,0.5-bh); // +/-0.5 not to hide chart border
							}
						}
						bar.attr({fill:catElem.color,stroke:catElem.color,opacity:CP.datasets[dsIdx].params.opacity})
						.data('js',catElem)
						.hover(function(e){setBarEmphasis(this,'on',e)},function(){setBarEmphasis(this,'off')})
						.click(function() {barElemClicked(this.data('js'))});
						catElem.svg=bar;
					}
					
					/** Point **/
					else if (CP.datasets[dsIdx].params.type==='point') {
						let pv=Math.round((toScale[axis](catElem.value)-minValue[axis])/pix2valScale[axis]),
						    px,py;
						if (CP.horizontalCP) {
							px=catPlotX0+pv+baseShift;
							py=groupCenterPos;
						}
						else { // default
							px=groupCenterPos;
							py=catPlotStartY0-pv-baseShift;
						}
						/* Point size calculation */					
						var pSize;
						if (typeof(CP.datasets[dsIdx].params.pointSizeRule)==='object') {
							let rule=CP.datasets[dsIdx].params.pointSizeRule,
							transformFunc=(rule.transform)? rule.transform.func : (x)=>x, // arrow function
							transR0=transformFunc(rule.range[0]),
							transR1=transformFunc(rule.range[1]),
							transRange=(transR0 <= transR1)? [transR0,transR1] : [transR1,transR0],
							mappedValue=(catElem.propValues && idvis.isNumber(catElem.propValues[rule.mapPropIndex]))? catElem.propValues[rule.mapPropIndex] : rule.range[0];
							mappedValue=Math.min(transRange[1],Math.max(transformFunc(mappedValue),transRange[0])); // make sure value is within range
							//pSize = Math.round( 0.5 * (groupSize * 0.1 + 0.9 * groupSize * (mappedValue-rule.range[0])/(rule.range[1]-rule.range[0])) ); // 0.5 <- radius
							pSize = Math.round(Math.max(3, 0.5 * (maxPointDiameter * (Math.sqrt(mappedValue)-Math.sqrt(transRange[0]))/(Math.sqrt(transRange[1])-Math.sqrt(transRange[0]))))); // 0.5 <- radius
						}
						else {pSize=CP.datasets[dsIdx].params.pointSize;}
//console.log('r',pSize);

						catElem.svg=idvis.lib.drawChartPoint(CP.canvas,px,py,pSize,CP.datasets[dsIdx].params.pointShape).attr({fill:catElem.color,stroke:'none',opacity:CP.datasets[dsIdx].params.opacity}) //,stroke:catElem.color
						.data('js',catElem)
						//.hover(function(){setPointEmphasis(this,'on')},function(){setPointEmphasis(this,'off')})
						//.click(function() {pointElemClicked(this.data('js'))});
						.hover(function(e){let type=(idvis.modKeyPressed || CP.datasets[dsIdx].params.type==='point')? 'max' : 'min'; idvis.lib.setLabelDisplay(this,'on',type,e);},function(){idvis.lib.setLabelDisplay(this,'off');})
						.click(function(){idvis.lib.setPointSelection(this,'auto','max');});
					}
				}
			}
		}
		
		/***** Point connection *****/
		for (let dsIdx=0; dsIdx<CP.datasets.length; dsIdx++) {
			if (!CP.datasets[dsIdx].dataConnection) continue;
			let pathPoints=[];
			for (let catIdx=0; catIdx<CP.categories.length; catIdx++) {
				if (!CP.datasets[dsIdx].data[catIdx]) {
					if (!CP.datasets[dsIdx].dataConnection.ignoreMissing) pathPoints.push(null); // Gap in path
					continue;
				}
				let svgObj=CP.datasets[dsIdx].data[catIdx].svg;
				if (CP.datasets[dsIdx].params.type=='point') {pathPoints.push(svgObj);}
				else {
					let x,y;
					if (CP.datasets[dsIdx].params.type=='box') {
						if (CP.horizontalCP) {
							x=svgObj[1].attr('path')[0][1]; // median
							y=svgObj[2].attr('path')[2][2]; // wiskers
						}
						else {
							x=svgObj[2].attr('path')[2][1]; // wiskers
							y=svgObj[1].attr('path')[0][2]; // median
						}
					}
					else { // bar/histogram
						if (CP.horizontalCP) {
							x=(svgObj.data('js').value >= 0)? svgObj.attr('x')+svgObj.attr('width') : svgObj.attr('x');
							y=svgObj.attr('y')+(svgObj.attr('height')/2);
						}
						else {
							x=svgObj.attr('x')+(svgObj.attr('width')/2);
							y=(svgObj.data('js').value >= 0)? svgObj.attr('y') : svgObj.attr('y')+svgObj.attr('height');
						}
					}
					pathPoints.push(CP.canvas.circle(x,y,10));
				}
			}
			if (pathPoints.length) {
				let connectPath=idvis.lib.drawPath(CP,CP.datasets[dsIdx],pathPoints); // chartLibrary
				if (CP.datasets[dsIdx].params.type === 'point') {
					let catElem;
					for (let i=0; i<CP.categories.length; i++) {
						if (CP.datasets[dsIdx].data[i]) {
							catElem=CP.datasets[dsIdx].data[i];
							break;
						}
					}
					connectPath.insertBefore(catElem.svg); // move behind first point
				}
				else {
					for (let i=0; i<pathPoints.length; i++) {pathPoints[i].remove();}
				}
				connectPath.data('js',CP.datasets[dsIdx])
				.hover(function(){this.attr({'stroke-width':this.attr('stroke-width')*2})},function(){this.attr({'stroke-width':this.data('js').dataConnection.width})});
				CP.datasets[dsIdx].dataConnection.svg=connectPath;
			}
		}

		/***** Features *****/
		if (CP.features) {
			for (let i=0; i<CP.features.length; i++) {
				let f=CP.features[i],
					anchorPoints=[];
//console.log(i+': '+f.type);
				if (CP.horizontalCP) {
					/* Threshold line */
					if (f.type=='line') {
						if (f.axis.match('X')) { // X or X2: Numerical axi(e)s
							f.pathPos=catPlotX0+Math.round((toScale[f.axis](f.value)-minValue[f.axis])/pix2valScale[f.axis]);
							//if (f.pathPos >= chartX && f.pathPos <= chartX0+chartW) { // if/when zoom implemented
								anchorPoints=[[f.pathPos,catPlotY0],[f.pathPos,catPlotY0+catPlotH]];
							//}
						}
						else { // Y: Category axis
							if (histoClass) { //CP.datasets[0].params.type==='histogram'
								f.pathPos=catPlotY0+Math.round((f.value-minValue[f.axis])/pix2valScale[f.axis]); // catPlotY0 (NOT catPlotStartY0) because "0" is at top of chart	
								//if (f.pathPos >= chartY && f.pathPos <= chartY0+chartH) { // if/when zoom implemented
									anchorPoints=[[catPlotX0,f.pathPos],[catPlotX0+catPlotW,f.pathPos]];
								//}
							}
							else { // value is a category index
								//TODO: draw "AFTER" the catory
							}
						}
					}
				}
				else { // Vertical (default)
					/* Threshold line */
					if (f.type==='line') {
						if (f.axis.match('Y')) { // Y or Y2: Numerical axi(e)s
							f.pathPos=catPlotStartY0-Math.round((toScale[f.axis](f.value)-minValue[f.axis])/pix2valScale[f.axis]);	
							//if (f.pathPos >= chartY && f.pathPos <= chartY0+chartH) { // if/when zoom implemented
								anchorPoints=[[catPlotX0,f.pathPos],[catPlotX0+catPlotW,f.pathPos]];
							//}
						}
						else { // X: Category axis
							if (histoClass) { //CP.datasets[0].params.type==='histogram'
								f.pathPos=catPlotX0+Math.round((f.value-minValue[f.axis])/pix2valScale[f.axis]);
								//if (f.pathPos >= chartX && f.pathPos <= chartX0+chartW) { // if/when zoom implemented
									anchorPoints=[[f.pathPos,catPlotY0],[f.pathPos,catPlotY0+catPlotH]];
								//}
							}
							else { // value is a category index
								//TODO: draw "AFTER" the catory
							}
						}
					}
				}
				idvis.lib.drawFeature(CP,f,anchorPoints);
			}
		}

		/***** Drawing scales legends *****/
		if (this.showScales && typeof(this.scales)==='object') {
			var refPosY=catPlotY0+7, scaleWidth=18, scaleHeight=(catWidth*this.categories.length > 250)? 100 : 50, fontSize=14; //axisTitleFontSize;  refPosX=catPlotX0+catPlotW+20,
			/** Color scale **/
			if (this.scales.color && this.scales.color.length) {
				this.scales.color.forEach((rule) => {
					if (!rule.range) return; // declared wo range and never matched data
					let ruleLabel=(rule.transform)? rule.transform.label : rule.map,
						tLabel=CP.canvas.text(legendStartX,refPosY,ruleLabel+':').attr({'font-size':fontSize,'text-anchor':'start','font-weight':'bold',fill:'#000'}),
						tBox=tLabel.getBBox(),
						startY0=refPosY + tBox.height/2 + 5,
						//gradient='90-rgb(0%,0%,100%)-rgb(70%,70%,0%)';
						gradient='90-'+scaleGradients[rule.gradient](0,1)+'-'+scaleGradients[rule.gradient](1,0);
					/* Scale */
					CP.canvas.rect(legendStartX,startY0,scaleWidth,scaleHeight).attr({stroke:'none',fill:gradient});
					/* Values and ticks */
					let transformFunc=(rule.transform)? rule.transform.func : (x)=>x, // arrow function
						transR0=transformFunc(rule.range[0]),
						transR1=transformFunc(rule.range[1]),
						transRange=(transR0 <= transR1)? [transR0,transR1] : [transR1,transR0],
						valueTicks=[transRange[0]],
						dividers=[4,2,4/3,1]
						startY=startY0+1; // To move 1st tick inside gradiant
					dividers.forEach(d => valueTicks.push(transRange[0]+(transRange[1]-transRange[0])/d));
					for (let i=0; i<valueTicks.length; i++)	{
						CP.canvas.path('M'+legendStartX+' '+startY+' l5 0 m'+(scaleWidth-10)+' 0 l5 0').attr({stroke:'#FFF'});
						CP.canvas.text(legendStartX+scaleWidth+2,startY,idvis.formatNumber(valueTicks[i])).attr({'font-size':9,'text-anchor':'start',fill:'#000'});
						startY=(i < valueTicks.length-2)? startY0 + scaleHeight/dividers[i] : startY0 + scaleHeight - 1;// To move last tick inside gradiant
					}
					refPosY+=scaleHeight+30;
				});
			}
			/* Point size scale */
			if (this.scales.pointSize && this.scales.pointSize.length) {
				let ruleIdx=-1;
				this.scales.pointSize.forEach((rule) => {
					if (!rule.range) return; // declared wo range and never matched data
					ruleIdx++;
					let ruleLabel=(rule.transform)? rule.transform.label : rule.map,
						tLabel=CP.canvas.text(legendStartX,refPosY,ruleLabel+':').attr({'font-size':fontSize,'text-anchor':'start','font-weight':'bold',fill:'#000'}),
						tBox=tLabel.getBBox(),
						transformFunc=(rule.transform)? rule.transform.func : (x)=>x, // arrow function
						transR0=transformFunc(rule.range[0]),
						transR1=transformFunc(rule.range[1]),
						transRange=(transR0 <= transR1)? [transR0,transR1] : [transR1,transR0],
						delta=transRange[1]-transRange[0],
						deltaRad=Math.sqrt(transRange[1])-Math.sqrt(transRange[0]),
						//valueRange=[Math.round((rule.range[1]-rule.range[0])/8),Math.round((rule.range[1]-rule.range[0])/4),Math.round((rule.range[1]-rule.range[0])/2),rule.range[1]],
						valueRange=[delta/8+transRange[0],delta/4+transRange[0],delta/2+transRange[0],transRange[1]],
						radiusRange=[];
					valueRange.forEach(v => radiusRange.push(Math.round(Math.max(3,0.5 * maxPointDiameter * (Math.sqrt(v)-Math.sqrt(transRange[0]))/deltaRad))));
					let startY=refPosY + tBox.height/2 + 5 + radiusRange[0],
						startX=legendStartX+radiusRange[valueRange.length-1];
					for (let i=0; i<valueRange.length; i++)	{
						// CP.canvas.circle(startX,startY,radiusRange[i]).attr({fill:categoryLabelColor || idvis.colorList[ruleIdx % idvis.colorList.length],stroke:'none'});
						idvis.lib.drawChartPoint(CP.canvas,startX,startY,radiusRange[i],scalesPointShapes[rule.map]).attr({fill:categoryLabelColor || idvis.colorList[ruleIdx % idvis.colorList.length],stroke:'none'});
						CP.canvas.text(startX+radiusRange[valueRange.length-1]+2,startY,idvis.formatNumber(valueRange[i])).attr({'font-size':9,'text-anchor':'start',fill:'#000'});
						startY+=(i < valueRange.length-1)? radiusRange[i]+radiusRange[i+1]+3 : radiusRange[i]+30;
					}
					refPosY+=scaleHeight+30;
				});
			}
		}

/* OBSOLETE -> this.drawDatasetIcon
		function drawIcon(type,x,y,color) { // x,y: top left corner of icon
			let svgSet=CP.canvas.set();
			if (type=='bar') {
				svgSet.push(CP.canvas.path('M'+x+' '+(y+10)+',l18 0'),
							CP.canvas.rect(x+2,y+4,4,6),
							CP.canvas.rect(x+7,y+10,4,4),
							CP.canvas.rect(x+12,y+2,4,8)
						   );
			}
			else if (type=='box') {
				svgSet.push(CP.canvas.rect(x+3,y+5,12,7),
							CP.canvas.path('M'+(x+3)+' '+(y+8)+' l12 0').attr({'stroke-width':2}),
							CP.canvas.path('M'+(x+5)+' '+(y+1)+' l8 0'),
							CP.canvas.path('M'+(x+9)+' '+(y+1)+' l0 4'),
							CP.canvas.path('M'+(x+5)+' '+(y+16)+' l8 0'),
							CP.canvas.path('M'+(x+9)+' '+(y+12)+' l0 4')
							);
			}
			else if (type=='point') {
				svgSet.push(CP.canvas.circle(x+9,y+7,6)); //.attr({stroke:'none',fill:color}));
			}
			else if (type=='histogram') {
				svgSet.push(CP.canvas.path('M'+x+' '+(y+13)+',l18 0'),
							CP.canvas.rect(x+2,y+9,4,4),
							CP.canvas.rect(x+7,y,4,13),
							CP.canvas.rect(x+12,y+5,4,8)
						   );
			}
			svgSet.attr({stroke:color,fill:color,'fill-opacity':0.5});
			return svgSet;
		}
*/
	};

	this.drawDatasetIcon = function(canvas,dataset) { // x,y: top left corner of icon
		let type=dataset.params.type,
			color=dataset.params.color,
			svgSet=canvas.set(),
			x=1,y=1; // adjustment variables
		if (type==='bar') {
			svgSet.push(canvas.path('M'+x+' '+(y+10)+',l18 0'),
						canvas.rect(x+2,y+4,4,6),
						canvas.rect(x+7,y+10,4,4),
						canvas.rect(x+12,y+2,4,8)
					   );
		}
		else if (type==='box') {
			svgSet.push(canvas.rect(x+3,y+5,12,7),
						canvas.path('M'+(x+3)+' '+(y+8)+' l12 0').attr({'stroke-width':2}),
						canvas.path('M'+(x+5)+' '+(y+1)+' l8 0'),
						canvas.path('M'+(x+9)+' '+(y+1)+' l0 4'),
						canvas.path('M'+(x+5)+' '+(y+16)+' l8 0'),
						canvas.path('M'+(x+9)+' '+(y+12)+' l0 4')
						);
		}
		else if (type==='point') {
			svgSet.push(canvas.circle(x+9,y+7,6)); //.attr({stroke:'none',fill:color}));
		}
		else if (type==='histogram') {
			svgSet.push(canvas.path('M'+x+' '+(y+13)+',l18 0'),
						canvas.rect(x+2,y+9,4,4),
						canvas.rect(x+7,y,4,13),
						canvas.rect(x+12,y+5,4,8)
					   );
		}
		svgSet.attr({stroke:color,fill:color,'fill-opacity':0.5});
		//return svgSet;
	}
	
	/*********************** Events management ******************************/
/*
	this.toggleDatasetVisibility = function(chkbox) { // public
		if (chkbox.readOnly) chkbox.checked=chkbox.readOnly=false;
		else if (!chkbox.checked) chkbox.readOnly=chkbox.indeterminate=true;
		let visStatus=(chkbox.checked)? 2 : (chkbox.indeterminate)? 1 : 0,
			dsIdx=chkbox.value;
		/ *Connection line* /
		let dsLine=CP.datasets[dsIdx].dataConnection;
		if (dsLine) {
			if (visStatus==0) {dsLine.svg.hide();}
			else if (visStatus==1) {dsLine.svg.attr({opacity:0.3,stroke:'#AAA'}).show();}
			else {dsLine.svg.attr({opacity:dsLine.opacity,stroke:dsLine.color}).show();}
		}
		/ *Points* /
		for (let catIdx=0; catIdx<CP.categories.length; catIdx++) {
			if (!CP.datasets[dsIdx].data[catIdx]) continue;
			let point=CP.datasets[dsIdx].data[catIdx];
			if (visStatus==0) {point.svg.hide();}
			else if (visStatus==1) {point.svg.attr({opacity:0.3,fill:'#AAA',stroke:'#AAA'}).show();}
			else {point.svg.attr({opacity:CP.datasets[dsIdx].params.opacity,fill:CP.datasets[dsIdx].params.color,stroke:CP.datasets[dsIdx].params.color}).show();}
			if (visStatus < 2 && point.svg.data('labelObj')) point.svg.data('labelObj').hide(); // unselected point
		}
	};
*/

	function setCatEmphasis(svgObj,action,e) {
		let catObj=svgObj.data('js');
		//let color=(action=='on')? '#F00' : catObj.color;
		if (action==='on') {
			svgObj.attr({fill:'#F00'});
			if (categoryOnClick) svgObj.attr({cursor:'pointer'});
			if (CP.horizontalCP) {highlightZone.attr({y:catPlotY0+(catWidth*catObj.index)});} else {highlightZone.attr({x:catPlotX0+(catWidth*catObj.index)});}
			highlightZone.show();
			//if (catObj.popupSvg) {catObj.popupSvg.show();}
			//else {
				let x=svgObj.attr('x'),
				    y=svgObj.attr('y');
				catObj.popupSvg=idvis.lib.drawLabel(CP.canvas,svgObj,x,y,1,1,catObj.popupText,catObj.labelColor,null,e);
			//}
		}
		else {
			svgObj.attr({fill:catObj.labelColor});
			if (categoryOnClick) svgObj.attr({cursor:'default'});
			catObj.popupSvg.remove(); //.hide();
			catObj.popupSvg=null;
			highlightZone.hide();
		}
	}
	function setBoxEmphasis(svgObj,action,e) {
		var boxObj=svgObj.data('js');
//let color=(action=='on')? '#F00' : boxObj.color;
		//boxObj.labelSvg.attr({fill:color});
		//boxObj.svg.attr({fill:color,stroke:color});
		//for (var i=0; i< boxObj.outlierList.length; i++) {
		//	if (!boxObj.outlierList[i].isSelected) {boxObj.outlierList[i].svg.attr({fill:color,stroke:color});}
		//}
		if (action==='on') {
			if (boxOnClick) svgObj.attr({cursor:'pointer'});
			if (CP.horizontalCP) {highlightZone.attr({y:catPlotY0+(catWidth*boxObj.catIndex)});} else {highlightZone.attr({x:catPlotX0+(catWidth*boxObj.catIndex)});}
			highlightZone.show();
			//if (boxObj.popupSvg) {boxObj.popupSvg.show();}
			//else {
				//let setStrg=(CP.datasets.length > 1)? '\n'+boxObj.dataset.params.name : '',
				//    text=boxObj.label+setStrg+':\n'+boxObj.pointList.length+' values\nUpper W='+boxObj.upperWhisker+'\nUpper Q='+boxObj.upperQuartile+'\nMedian='+boxObj.median+'\nLower Q='+boxObj.lowerQuartile+'\nLower W='+boxObj.lowerWhisker+'\nMean='+(boxObj.mean.toPrecision(4)*1)+'\n'+boxObj.numOutliers+' outlier(s)';
				let text=boxObj.info();
				boxObj.popupSvg=idvis.lib.drawLabel(CP.canvas,svgObj,boxObj.x,boxObj.y,1,1,text,boxObj.color,null,e);
			//}
		}
		else {
			if (boxOnClick) svgObj.attr({cursor:'default'});
			boxObj.popupSvg.remove(); // hide();
			boxObj.popupSvg=null;
			highlightZone.hide();
		}
	}
	function setBarEmphasis(svgObj,action,e) {
		barObj=svgObj.data('js');
		//var color=(action=='on')? '#F00' : barObj.color;
		if (action==='on') {
			if (barOnClick) svgObj.attr({cursor:'pointer'});
			let x,y;
			if (CP.horizontalCP) {
				highlightZone.attr({y:catPlotY0+spaceBeforeBins+(catWidth*barObj.catIndex)});
				x=(barObj.value >= 0)? svgObj.attr('x')+svgObj.attr('width') : svgObj.attr('x');
				y=svgObj.attr('y')+(svgObj.attr('height')/2);
			}
			else { // default
				highlightZone.attr({x:catPlotX0+spaceBeforeBins+(catWidth*barObj.catIndex)});
				x=svgObj.attr('x')+(svgObj.attr('width')/2);
				y=(barObj.value >= 0)? svgObj.attr('y') : svgObj.attr('y')+svgObj.attr('height');
			}
			highlightZone.show();
			//if (barObj.popupSvg) {barObj.popupSvg.show();}
			//else {
				let setStrg=(CP.datasets.length > 1)? '\n'+barObj.dataset.params.name : '',
					valueLabel=CP.axesTitles[barObj.dataset.params.axis].shortTitle || 'Value';
				barObj.popupSvg=idvis.lib.drawLabel(CP.canvas,svgObj,x,y,1,1,barObj.label+setStrg+'\n'+valueLabel+': '+barObj.value,barObj.color,null,e);
			//}
		}
		else {
			if (barOnClick) svgObj.attr({cursor:'default'});
			barObj.popupSvg.remove(); //.hide();
			barObj.popupSvg=null;
			highlightZone.hide();
		}
	}
/* Handled by idvis
	function setPointEmphasis(svgObj,action) {
		pointObj=svgObj.data('js');
		//var color=(action=='on')? '#F00' : pointObj.color;
		if (action=='on') {
			if (CP.horizontalCP) {
				highlightZone.attr({y:catPlotY0+(catWidth*pointObj.catIndex)});
			}
			else { // default
				highlightZone.attr({x:catPlotX0+(catWidth*pointObj.catIndex)});
			}
			highlightZone.show();
			if (pointObj.popupSvg) {pointObj.popupSvg.show();}
			else {
				let setStrg=(CP.datasets.length > 1)? '\n'+pointObj.dataset.params.name : '';
				pointObj.popupSvg=idvis.lib.drawLabel(CP.canvas,svgObj,svgObj.attr('cx'),svgObj.attr('cy'),svgObj.attr('r')+1,svgObj.attr('r')+1,pointObj.label+setStrg+'\nValue: '+pointObj.value,pointObj.color);
			}
		}
		else {
			pointObj.popupSvg.hide();
			highlightZone.hide();
		}
	}
*/
	function boxElemClicked(boxObj) {
		if (idvis.modKeyPressed && excludeOnClick) {
			if (excludeOnClick(boxObj.externalID,boxObj.excluded)) { // 2 parameters boxID,isExcluded
				let update=(boxObj.excluded)? [false,'',''] : [true,'- ','. '];
				boxObj.excluded=update[0];
				boxObj.svg.attr({'stroke-dasharray':update[1]});
				for (let j=0; j<boxObj.pointList.length; j++) {
					boxObj.pointList[j].svg.attr({'stroke-dasharray':update[2]});
				}
			}
			else {alert('ERROR: External resource indicates update failure!');}
		}
		else if (boxOnClick) {boxOnClick(boxObj.dataset.params.externalID,boxObj.label,boxObj.externalID,CP.categories[boxObj.catIndex].externalID);}
	}

	function barElemClicked(barObj) {
		if (barOnClick) {barOnClick(barObj.dataset.params.externalID,barObj.label,barObj.externalID,CP.categories[barObj.catIndex].externalID);}
	}
/*	
	function pointElemClicked(pointObj) {
		if (pointOnClick) {pointOnClick(pointObj.externalID);}
	}
*/	
/*
	this.changeThreshold=function(newValue) {
		var scaledNewValue=toScale[axis](newValue);
		if (!idvis.isNumber(scaledNewValue)) { // chartLibrary
			alert("'"+newValue+"' is not a valid number!");
			document.getElementById(thresholdInputID).value=thresholdValue;
			return;
		}
		if (scaledNewValue < minValue[axis] || scaledNewValue > maxValue) {
			alert("'"+newValue+"' is out of range!");
			document.getElementById(thresholdInputID).value=thresholdValue;
			return;
		}
		cpThresholdLine.setValues(newValue);
		thresholdValue=newValue;
		//var shift=Math.round((toScale[axis](newValue)-toScale[axis](cpThresholdLine.startValue))/pix2valScale[axis]);
		if (cpThresholdLine.axis=='X') { // horizontal CP
			var shift=Math.round((scaledNewValue-toScale[axis](cpThresholdLine.startValue))/pix2valScale[axis]);
			cpThresholdLine.path.transform('t'+shift+',0');
		}
		else { // vertical CP
			var shift=Math.round((toScale[axis](cpThresholdLine.startValue)-scaledNewValue)/pix2valScale[axis]);
			cpThresholdLine.path.transform('t0,'+shift);
		}
	};

	this.getExcludedBoxes=function() {
		var excludedBoxes=[];
		for (var i=0; i<categories.length; i++) {
			if (categories[i].median < cpThresholdLine.initValue) excludedBoxes.push(categories[i].externalID);
		}
		if (!cpData.threshold.action[1](cpThresholdLine.initValue,excludedBoxes)) {
			alert('ERROR: External resource indicates update failure!');
			return;
		}
		// OK to exclude
		for (var i=0; i<categories.length; i++) {
			var update=[];
			if (categories[i].median < cpThresholdLine.initValue) {
				if (!categories[i].excluded) {
					update=[true,'- ',0.5,'. '];
				}
			}
			else if (categories[i].excluded) { // median >= threshold
				update=[false,'',1,''];
			}
			if (update.length) {
				categories[i].excluded=update[0];
				categories[i].svg.attr({'stroke-dasharray':update[1],'stroke-opacity':update[2]});
				for (var j=0; j<categories[i].pointList.length; j++) {
					categories[i].pointList[j].svg.attr({'stroke-dasharray':update[3],'stroke-opacity':update[2]});
				}
			}
		}
	};
*/

/*
	function setBoxPointEmphasis(boxPoint,action) {
		if (action=='on') {
			if (boxPoint.popupSvg) {boxPoint.popupSvg.show();}
			else {
				let d=boxPoint.svg.attr('r')+1,
				    text=(boxPoint.label)? boxPoint.label+': ' : '';
				text+=boxPoint.value;
				boxPoint.popupSvg=idvis.lib.drawLabel(CP.canvas,boxPoint.svg,boxPoint.svg.attr('cx'),boxPoint.svg.attr('cy'),d,d,text,boxPoint.catElement.color);
			}
			boxPoint.svg.attr({fill:'#F00',stroke:'#F00'});
		}
		else if (!boxPoint.isSelected) { // off
			boxPoint.popupSvg.hide();
			boxPoint.svg.attr({fill:boxPoint.catElement.color,stroke:boxPoint.catElement.color});
		}
	}

	function selectBoxPoint(boxPoint) { // on/off switch
		//if (action=='auto') { // invert selection
			if (boxPoint.isSelected) { // already showing label => set OFF
				boxPoint.isSelected=false;
				setBoxPointEmphasis(boxPoint,'off');
			}
			else { // set ON
				boxPoint.isSelected=true;
				setBoxPointEmphasis(boxPoint,'on');
			}
		//}
/-*
		else if (action=='on') { // set ON
			if (!outObj.isSelected) {
				outObj.isSelected=true;
				setBoxPointEmphasis(outObj,'on');
			}
		}
		else { // OFF
			outObj.isSelected=false;
			setBoxPointEmphasis(outObj,'off');
		}
*-/
	}
 */

    this.showBoxValues=function(status) {
        for (let d=0; d < CP.datasets.length; d++) {
            if (CP.datasets[d].params.type != 'box') {continue;}
            for (let c=0; c < CP.categories.length; c++) {
                if (!CP.datasets[d].data[c]) {continue;}
                for (let i=0; i < CP.datasets[d].data[c].pointList.length; i++) {
                    var boxPoint=CP.datasets[d].data[c].pointList[i];
                    if (boxPoint.isOutlier || boxPoint.isExcluded) {continue;} // outliers & excluded always visible
                    if (status===true) {boxPoint.svg.show();}
                    else {
                        //if (boxPoint.isSelected) {selectBoxPoint(boxPoint);} // on/off switch
						if (boxPoint.svg.data('showLabel')) {idvis.lib.setPointSelection(boxPoint.svg,'off','max');} // on/off switch
                        boxPoint.svg.hide();
                    }
                }
            }
        }
    };
	
	
	/*====================================== Nested objects ===================================*/
	
	/*****************************************************************
					  Box object
	******************************************************************/
	var box = function(dataset,catIdx,boxInfo,baseValue) {
		//this.dsIndex=dsIdx;
		this.dataset=dataset;
		this.catIndex=catIdx;
		this.label=boxInfo.label;
		this.externalID=boxInfo.id || catIdx;
		this.excluded=(boxInfo.excluded)? true : false;
		this.color=(boxInfo.color)? boxInfo.color : (this.dataset.params.chart.categories[catIdx].color)? this.dataset.params.chart.categories[catIdx].color : dataset.params.color;
		this.numValues=null;
		this.median=null;
		this.mean=null;
		this.lowerQuartile=null;
		this.upperQuartile=null;
		this.lowerWhisker=null;
		this.upperWhisker=null;
		this.extremeValues=[];
	//this.outlierList=[];
		this.numOutliers=0;
		this.numExcluded=0;
		this.pointList=[];
		this.computeBox(boxInfo);
		this.x=this.y=null;
		this.baseValue=baseValue;
		//this.labelSvg; // graphical label
		this.svg=null; // set of graphical objects
		this.meanSvg=null;
		this.popupSvg=null;
	};
	box.prototype = {
		computeBox: function(boxInfo) {
			let tmp2dArray=[],
				excludArr=[],
			    numErrors=0;
			for (let i=0; i<boxInfo.data.length; i++) { // converting into a 2D array
				let arr=(boxInfo.data[i]+'').split(':'); // to string (fails if true number)
				arr[0]*=1;
				if (idvis.isNumber(arr[0])) {
					if (arr[2]) {excludArr.push(arr)} // exclusion flag
					else {tmp2dArray.push(arr);}
				}
				else {numErrors++;}
			}
			if (numErrors>0) {alert('WARNING: '+numErrors+' non-valid value(s) found for '+this.label);}
			this.numValues=tmp2dArray.length;
			var dataList=tmp2dArray.sort(function(a,b){return a[0]-b[0];}); // ascending order
			this.extremeValues=[dataList[0][0],dataList[dataList.length-1][0]];
			//Box
			var medianRk=this.numValues/2, minMedianRk=Math.floor(medianRk), maxMedianRk=Math.ceil(medianRk);
			this.median=(medianRk != minMedianRk)? dataList[minMedianRk][0] : (dataList[medianRk-1][0]+dataList[medianRk][0])/2;
			var boxMinRk=maxMedianRk/2, minBoxMinRk=Math.floor(boxMinRk);
			this.lowerQuartile=(boxMinRk != minBoxMinRk)? dataList[minBoxMinRk][0] : (dataList[boxMinRk-1][0]+dataList[boxMinRk][0])/2;
			var boxMaxRk=this.numValues-boxMinRk, minBoxMaxRk=Math.floor(boxMaxRk);
			this.upperQuartile=(boxMaxRk != minBoxMaxRk)? dataList[minBoxMaxRk][0] : (dataList[boxMaxRk-1][0]+dataList[boxMaxRk][0])/2;
			this.mean=0;
			for (let i=0; i<this.numValues; i++) {this.mean+=dataList[i][0];}
			this.mean/=this.numValues;
//console.log(dataList.join(' '));
//console.log(medianRk+':'+minMedianRk+', '+boxMinRk+':'+minBoxMinRk+', '+boxMaxRk+':'+minBoxMaxRk);
			//Whiskers
			if (this.dataset.params.chart.outlierRule=='IQR') {
				var maxWhkSize=(this.upperQuartile-this.lowerQuartile)*1.5;
				//Lower
				var minWhkValue=this.lowerQuartile-maxWhkSize;
				var maxWhkValue=this.upperQuartile+maxWhkSize;
				for (let i=0; i<this.numValues; i++) {
					let isOutlier=false;
					if (dataList[i][0] < minWhkValue) {isOutlier=true;} // lower outlier
					else {
						if (this.lowerWhisker===null) {this.lowerWhisker=dataList[i][0];}
						if (dataList[i][0] > maxWhkValue) {isOutlier=true;} // upper outlier
						else {this.upperWhisker=dataList[i][0];} // takes all values until dataList[i][0] > maxWhkValue
					}
					this.pointList.push(new boxPoint(this,dataList[i],isOutlier));
					if (isOutlier===true) this.numOutliers++;
				}
			}
			else { // eg. 0.09 or 0.02
				//Lower
				var firstWhkIdx=Math.ceil(this.numValues*CP.outlierRule)-1;
				this.lowerWhisker=dataList[firstWhkIdx][0];
				//upper
				var lastWhkIdx=this.numValues-firstWhkIdx-1;
//console.log(firstWhkIdx,lastWhkIdx);
				this.upperWhisker=dataList[lastWhkIdx][0];          
				for (let i=0; i<this.numValues; i++) {
					let isOutlier=(i < firstWhkIdx || i > lastWhkIdx)? true : false;
					this.pointList.push(new boxPoint(this,dataList[i],isOutlier,false));
					if (isOutlier===true) this.numOutliers++;
				} 
			}
			this.numExcluded=excludArr.length;
			for (let i=0; i<this.numExcluded; i++) {
				this.pointList.push(new boxPoint(this,excludArr[i],false,true));
			}
		},
		info: function() {
			let setStrg=(this.dataset.params.chart.datasets.length > 1)? '\n'+this.dataset.params.name : '',
				infoStrg=this.label+setStrg+':\n'
					+this.pointList.length+' values\nUpper W='+idvis.formatNumber(this.upperWhisker)
					+'\nUpper Q='+idvis.formatNumber(this.upperQuartile)
					+'\nMedian='+idvis.formatNumber(this.median)
					+'\nLower Q='+idvis.formatNumber(this.lowerQuartile)
					+'\nLower W='+idvis.formatNumber(this.lowerWhisker)
					+'\nMean='+idvis.formatNumber(this.mean)
					+'\n'+this.numOutliers+' outlier(s)';
			if (this.numExcluded) {
				infoStrg+='\n'+this.numExcluded+' excluded point(s)';
			}
			return infoStrg;		
		}
	};
	
	
	/*****************************************************************
					 Point object for box
	******************************************************************/
	var boxPoint = function(box,pointData,isOutlier,isExcluded=false) {
		this.catElement=box;
		this.catIndex=box.catIndex;
		this.dataset=box.dataset;
		this.value=pointData[0]*1;
		this.label=pointData[1] || null;
		this.externalID=pointData[2] || this.label;
		this.isOutlier=isOutlier;
		this.isExcluded=isExcluded;
		this.color=(isExcluded && boxPointExclusion.color)? boxPointExclusion.color : null,
		//this.x=this.y=null;
		this.svg=null; // svg object
		this.popupSvg=null;
		//this.isSelected=false;
	};
	boxPoint.prototype =  {
		info: function(type) {
			let infoStrg=this.label || this.dataset.params.chart.categories[this.catIndex].label || 'No label';
			if (type == 'min') {
				infoStrg+='\n'+idvis.formatNumber(this.value); //.toFixed(3);
			}
			else {
				infoStrg='\nBox: '+this.catElement.label;
				if (this.dataset.params.chart.datasets.length > 1) infoStrg+='\nSeries: '+this.dataset.params.name;
				let valueLabel=this.dataset.params.chart.axesTitles[this.dataset.params.axis].shortTitle || 'Value';
				infoStrg+='\n'+valueLabel+': '+idvis.formatNumber(this.value); //.toFixed(3);
			}
			if (this.isOutlier) {infoStrg+='\nOutlier';}
			else if (this.isExcluded) {infoStrg+=(boxPointExclusion.label)? '\n'+boxPointExclusion.label : '\nExcluded';}
			return infoStrg;
		}//,
/*
		getX: function() {
			return this.x;
		},
		getY: function() {
			return this.y;
		}
*/
	};
	
};


/* Random gaussian distribution */
idvis.ziggurat = function(){

  var jsr = 123456789;

  var wn = Array(128);
  var fn = Array(128);
  var kn = Array(128);

  function RNOR(){
    var hz = SHR3();
    var iz = hz & 127;
    return (Math.abs(hz) < kn[iz]) ? hz * wn[iz] : nfix(hz, iz);
  }

  this.nextGaussian = function(){
    return RNOR();
  };

  function nfix(hz, iz){
    var r = 3.442619855899;
    var r1 = 1.0 / r;
    var x;
    var y;
    while(true){
      x = hz * wn[iz];
      if( iz == 0 ){
        x = (-Math.log(UNI()) * r1); 
        y = -Math.log(UNI());
        while( y + y < x * x){
          x = (-Math.log(UNI()) * r1); 
          y = -Math.log(UNI());
        }
        return ( hz > 0 ) ? r+x : -r-x;
      }

      if( fn[iz] + UNI() * (fn[iz-1] - fn[iz]) < Math.exp(-0.5 * x * x) ){
         return x;
      }
      hz = SHR3();
      iz = hz & 127;
 
      if( Math.abs(hz) < kn[iz]){
        return (hz * wn[iz]);
      }
    }
  }

  function SHR3(){
    var jz = jsr;
    var jzr = jsr;
    jzr ^= (jzr << 13);
    jzr ^= (jzr >>> 17);
    jzr ^= (jzr << 5);
    jsr = jzr;
    return (jz+jzr) | 0;
  }

  function UNI(){
    return 0.5 * (1 + SHR3() / -Math.pow(2,31));
  }

  function zigset(){
    // seed generator based on current time
    jsr ^= new Date().getTime();

    var m1 = 2147483648.0;
    var dn = 3.442619855899;
    var tn = dn;
    var vn = 9.91256303526217e-3;
    
    var q = vn / Math.exp(-0.5 * dn * dn);
    kn[0] = Math.floor((dn/q)*m1);
    kn[1] = 0;

    wn[0] = q / m1;
    wn[127] = dn / m1;

    fn[0] = 1.0;
    fn[127] = Math.exp(-0.5 * dn * dn);

    for(var i = 126; i >= 1; i--){
      dn = Math.sqrt(-2.0 * Math.log( vn / dn + Math.exp( -0.5 * dn * dn)));
      kn[i+1] = Math.floor((dn/tn)*m1);
      tn = dn;
      fn[i] = Math.exp(-0.5 * dn * dn);
      wn[i] = dn / m1;
    }
  }
  zigset();
};
/*
To use it, just instantiate an instance of Ziggurat(), and call nextGaussian() on it to return a normally distributed random number with mean 0 and standard deviation 1:

var z = new Ziggurat();
z.nextGaussian();

*/





/*
####>Revision history<####
# 1.2.0 [FEATURE] Support for point shapes & improved popup display (PP 06/05/21)
# 1.1.0 [FEATURE] Added support for color/point size scales, properties and click event on category label (PP 03/12/20)
# 1.0.0 First stable version (PP 01/09/20)
# 0.0.1 Started from boxPlot 1.0.5 (PP 23/06/15)

*/