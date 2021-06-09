/*
################################################################################
# idvis.js         2.2.0                                                       #
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



/***************** Mouse and Keyboard events ********************/
document.onkeydown=function(e){
	let kbEv=window.event || e;
	if (kbEv.keyCode==16 || kbEv.keyCode==17 || kbEv.keyCode==18) { // ctrl or alt
		idvis.modKeyPressed=true;
	}
};
document.onkeyup=function(e){
	let kbEv=window.event || e;
	if (kbEv.keyCode==16 || kbEv.keyCode==17 || kbEv.keyCode==18) { // ctrl or alt
		idvis.modKeyPressed=false;
	}
};

/*****************************************************
     Top object with public attributes & methods
******************************************************/
var idvis = {
	author : 'Patrick Poullet',
	version : '2.2.0',
	license : 'CeCILL',

	registeredCharts : [], // list of chart drawn
	colorList : ['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B'],
	highlightColor: '#F00',	
	dimColor: '#AAA',
	dimSetOpacity: 0.5, // for dimming entire dataset
	dimPointOpacity: 0.2, // for dimming points based on thresholds
	pathHoverWidth: 4,
	modKeyPressed : false,
	closePopupIcon : null,

    isNumber : function(n) {
      return !isNaN(parseFloat(n)) && isFinite(n);
    },
	formatNumber : function(n,d) {
		if (!idvis.isNumber(n)) return n;
		n*=1;
		let n1=Math.abs(n),
		    //n1=(n0 > 1)? n0 : 1/n0;
			hasFixed=(idvis.isNumber(d))? true : false,
			d1=(hasFixed)? d : 2, // d could be '0'
		    nf=((n1 != 0 && n1 < 10**(-d1+1)) || n1 > 1000)? n.toExponential(d1) : (hasFixed || n1.toString().length > 4)? 1*(n.toFixed(d1)) : n;
/*		
		if (n1 >= 1000 || n0 < 0.1) {nf=n.toExponential(d);}
		else if (n1 > 1) {nf=n.toExponential(d);}
		if (n1 > 0.1 && n1 < 100) {nf=n.fixed(d);}
		if (n1.toString().length > 4) { // large number of many decimal digits
			
			nf=(n < 0.01 || n > 100 )? n.toExponential(d) : n.toFixed(2);
		}
*/
		return nf;
	},
	valuesProperties : function(values) {
		values[0]*=1;
		let prop={mean:values[0], min:values[0], max:values[0]};
		for (let i=1; i<values.length; i++) {
			values[i]*=1;
			prop.mean+=values[i];
			prop.min=Math.min(prop.min,values[i]);
			prop.max=Math.max(prop.max,values[i]);
		}
		prop.mean/=values.length;
		return prop;
	},
	sortNumber : function(a,b) {return a-b;},
	//sortNumberAsc : sortNumber,
	sortNumberDesc : function(a,b) {return b-a;},

	/************* Mouse position *************/
	getMousePositionInElement : function(e) {
		let target = e.target || e.srcElement,
			style = target.currentStyle || window.getComputedStyle(target, null),
			borderLeftWidth = parseInt(style.borderLeftWidth, 10),
			borderTopWidth = parseInt(style.borderTopWidth, 10),
			rect = target.getBoundingClientRect();
		return [e.clientX - borderLeftWidth - rect.left,e.clientY - borderTopWidth - rect.top];
	},

	/************* String processing *************/
	shortenText : function(myText,maxLength) {
		let textLength=myText.length;
		if (textLength <= maxLength) return myText;
		let subTextSize=Math.max(2,Math.floor((maxLength-3)/2));
		return myText.substring(0,subTextSize)+'...'+myText.substring(textLength-subTextSize);
	},
	ucFirst : function(myText) {
		myText+='';
		return myText.charAt(0).toUpperCase() + myText.slice(1);
	}

};

/*================================ idvis generic objects =============================*/

/************************ feature object ********************/
idvis.feature = function(C,fea) {
	this.isFeature=true;
	this.chart=C;
	this.name=fea.label; //used if editable
	this.color=fea.color || '#000';
	this.keepAbove1=(fea.keepAbove1)? true : false;
	this.pattern=fea.pattern || ''; // null for text
	this.width=fea.width || 1; // null for text
	this.opacity=fea.opacity || 1;
	this.type=fea.type; // line,function,custom,(path),text
	//if (!fea.axis) {this.axis='X';}
	//else if (fea.type !== 'line') {this.axis=(fea.axis.match('2'))? 'X2' : 'X';}
	//else {this.axis=fea.axis;}
	if (fea.type==='line') {this.axis=fea.axis || 'X';} // No axisY for 'line'
	else {
		this.axis=((fea.axis && fea.axis.match('2')) || (fea.axisX && fea.axisX.match('2')))? 'X2' : 'X';
		this.axisX=this.axis;
		this.axisY=(fea.axisY && fea.axisY.match('2'))? 'Y2' : 'Y';
	}
	this.path=null;
	this.popup=null;
	if (fea.type==='line') {
		this.lineType='line';
		this.startValue=fea.value;
		//this.roundValue=(fea.roundValue)? true : false; //optional (rounding level)
		this.editable=fea.editable; //optional
		this.baseLabel=fea.label;
		this.setValues(fea.value);
		this.pathPos=null; // needed when editable
	}
	else if (fea.type==='text') {
		this.label=fea.label;
		this.editable=false;
		this.x=fea.value[0] || 0;
		this.y=fea.value[1] || 0;
		this.anchor=fea.anchor || 'middle';
		this.size=fea.size || 12; // text size in px
		this.weight=fea.weight || 'normal';
	}
	else {
		this.label=fea.label;
		this.editable=false;
		this.lineType=fea.lineType || 'line';
		this.data=fea.value; // function if type='function'
		if (this.type==='custom') { // data is array of [x,y] arrays
			// find min/max x/y
			this.minX=this.maxX=this.minY=this.maxY=null;
			for (let i=0; i<this.data.length; i++) {
				if (this.minX===null) {
					this.minX=this.maxX=this.data[i][0];
					this.minY=this.maxY=this.data[i][1];
				}
				else {
					this.minX=Math.min(this.minX,this.data[i][0]);
					this.maxX=Math.max(this.maxX,this.data[i][0]);
					this.minY=Math.min(this.minY,this.data[i][1]);
					this.maxY=Math.max(this.maxY,this.data[i][1]);
				}
			}
		}
		else if (fea.range) { // function with range
			this.minX=fea.range[0];
			this.maxX=fea.range[1];
			if (fea.range.length==4) {
				this.minY=fea.range[2];
				this.maxY=fea.range[3];
			}
			else {
				this.minY=this.maxY=null;
				let dx=(this.maxX-this.minX)/20,
					x=this.minX;
				while (x < this.maxX+dx/2) { // to compensate for JavaScript poor number rounding
					if (this.minY===null) {
						this.minY=this.maxY=this.data(x);
					}
					else {
						let y=this.data(x);
						this.minY=Math.min(this.minY,y);
						this.maxY=Math.max(this.maxY,y);
					}
					x+=dx;
				}
			}
		}
	}
};
idvis.feature.prototype = { // prototype because can be called again if user changes feature's value (only for type='line')
	setValues: function(initValue) {
		this.initValue=initValue;
		this.value=(this.chart.convertValue)? this.chart.convertValue(this.axis,initValue) : initValue;
		if (this.baseLabel) {
			let dispValue;
			if (initValue < 1 && this.keepAbove1) {
				//dispValue=(this.roundValue)? '1/'+(1/initValue).toFixed(this.roundValue) : '1/'+(1/initValue);
				dispValue='1/'+idvis.formatNumber(1/initValue);
			}
			else {
				//dispValue=(this.roundValue)? (initValue*1).toFixed(this.roundValue) : initValue;
				dispValue=idvis.formatNumber(initValue);
			}
			this.label=this.baseLabel+': '+dispValue;
		}
		else {this.label=null;}
	}
};

/*================================ Library of idvis functions to manipulate charts =============================*/

idvis.lib = (function() {

	/******************* Register chart into idvis *********************/
	var registerChart = function(Chart) { // *public*
		idvis.registeredCharts.push(Chart);
		if (Chart.div) {
			Chart.div.innerHTML='<DIV style="padding:20px;font-weight:bold">Preparing chart. Please wait...</DIV>';
		}
		return idvis.registeredCharts.length-1; // index of chart in array
	};

	/******************* Dataset creation *********************/
	var addDataset = function(mainC,dsIdx,set) { // scatter plots only!
		let [dsPointShape,dsColor]=generateDatasetAES(mainC,dsIdx,set);
		mainC.datasets[dsIdx] = {
			data : [],
			//line : null, // SVG path
			params : {
				chart : mainC,
				index : dsIdx,
				axisX : (set.axisX && set.axisX.match('2'))? 'X2' : 'X', // dual X axes
				axisY : (set.axisY && set.axisY.match('2'))? 'Y2' : 'Y', // dual Y axes
				name : set.name || 'Dataset '+(dsIdx+1),
				pointShape: dsPointShape, //_decodePointShape(set.pointShape),
				color : dsColor, //set.color || idvis.colorList[dsIdx % idvis.colorList.length],
				opacity : set.opacity || mainC.pointOpacity || 1,
				pointLabelX: (set.pointLabels && set.pointLabels.axisX)? set.pointLabels.axisX : 'X',
				pointLabelY: (set.pointLabels && set.pointLabels.axisY)? set.pointLabels.axisY : 'Y',
				pointSizeName : (set.pointLabels && set.pointLabels.size)? set.pointLabels.size : null,
				//pointSizeName : set.pointSizeName || null,
				pointSizeRule : set.pointSizeRule || {min:1,max:20,scale:2},
				visible : 2,
				//line : set.connection || null, moved to dataConnection
				multiValuePoint : set.multiValuePoint || {linked:false,size:3}
			},
			dataConnection : set.connection || set.line || null // JS object with {type,color,pattern,width,opacity}
		};
		if (!mainC.datasets[dsIdx].params.pointSizeName && set.pointSizeName) {mainC.datasets[dsIdx].params.pointSizeName=set.pointSizeName;} // back compatibility for dset using pointSizeName instead of pointLabels
		/* Checking mandatory params */
		if (mainC.datasets[dsIdx].dataConnection) {
			if (!mainC.datasets[dsIdx].dataConnection.color) {mainC.datasets[dsIdx].dataConnection.color=mainC.datasets[dsIdx].params.color;}
			if (!mainC.datasets[dsIdx].dataConnection.opacity) {mainC.datasets[dsIdx].dataConnection.opacity=mainC.datasets[dsIdx].params.opacity;}			
			if (!mainC.datasets[dsIdx].dataConnection.width) {mainC.datasets[dsIdx].dataConnection.width=1;}			
		}
		if (mainC.datasets[dsIdx].params.multiValuePoint.linked===undefined) {mainC.datasets[dsIdx].params.multiValuePoint.linked=false;}
		if (!mainC.datasets[dsIdx].params.multiValuePoint.size) {mainC.datasets[dsIdx].params.multiValuePoint.size=5;}	
	};

	/******************* HTML interface *********************/
	var drawChartLayout = function (mainC,formPosition,canvasW,canvasDivID,formDivID) { // *public*
		// formPosition: top, bottom, right (default), left or none (hidden) 
		formPosition=(formPosition && formPosition.toString().toLowerCase().match('^(t|b|r|l|n)'))? formPosition.toString().toLowerCase() : 'right';
		let layoutHTML;
		if (formPosition.match('^t')) {
			layoutHTML='<TABLE><TR><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD></TR></TABLE>';
		}
		else if (formPosition.match('^b')) {
			layoutHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD></TR><TR><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		}
		else if (formPosition.match('^l')) {
			layoutHTML='<TABLE><TR><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD></TR></TABLE>';
		}
		else { // right|none
			let styleStrg=(formPosition.match('^n'))? 'display:none' : 'min-width:200px'; // hidden form
			layoutHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="'+styleStrg+'"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
		}
		mainC.div.innerHTML=layoutHTML;
		mainC.formPosition=formPosition;
	};

	var initializeForm = function (mainC,customFormElements=[]) { // *public*
		//customFormElements:
		// [0,10,20...]: control idvis default form elements =[isDisplayed,custom label for default field]
		// [all other indexes]: =custom html code to be injected in form
		if (mainC.chartID===undefined) { // just to be safe
			mainC.chartID=registerChart(mainC);
		}

		// const formName='idvisChartForm_'+mainC.divID;
		// let htmlString='<FORM name="'+formName+'" style="margin:0px">';
		let htmlString='';
		for (let i=0; i<10; i++) { // custom elements [0..9]
			if (customFormElements[i]) {htmlString+=customFormElements[i];}
		}

		/* Datasets */
		var pointDatasets=[],
			datasetIcons=[]
			datasetsLabel=mainC.datasetsLabel || 'Datasets';
		if (mainC.datasets.length) { // >=1 dataset
			//if (mainC.datasets.length > 1) datasetsLabel+='s';

			/* Datasets visibility [index=10] */
			if (mainC.showDatasets || (customFormElements[10] && customFormElements[10][0]===true) || (!customFormElements[10] && mainC.datasets.length > 1)) {
				if (customFormElements[10] && customFormElements[10][1]) datasetsLabel=customFormElements[10][1];
				htmlString+='<FIELDSET style="padding:2px;white-space:nowrap;max-height:200px;overflow:auto"><LEGEND><B>'+datasetsLabel+':</B></LEGEND>';
				let dsHtmlString='',
					numChkBoxes=0;
				for (let i=0; i<mainC.datasets.length; i++) {
					if (i > 0) {dsHtmlString+='<BR>';}
					if (mainC.datasets[i].params.selectable === undefined || mainC.datasets[i].params.selectable === true) {
						dsHtmlString+='<INPUT type="checkbox" name="datasetCheckboxes_'+mainC.chartID+'" id="datasetChkbox_'+mainC.chartID+'_'+i+'" value="'+i+'" onclick="idvis.lib.setDatasetDisplay(idvis.registeredCharts['+mainC.chartID+'],this)" checked>';
						numChkBoxes++;
					}
					// if (!mainC.datasets[i].params.type || (mainC.datasets[i].params.type === 'point' && !mainC.datasets[i].params.colorRule)) {pointDatasets.push(i);}
					// else {
						if (!mainC.datasets[i].params.type || (mainC.datasets[i].params.type === 'point' && !mainC.datasets[i].params.colorRule)) {pointDatasets.push(i);}
						const iconSpId=mainC.formDiv.id+'_ds'+i;
						dsHtmlString+='<SPAN id="'+iconSpId+'" style="vertical-align:top"></SPAN>';
						datasetIcons.push([i,iconSpId]);
					// }
					dsHtmlString+='<FONT id="dataset_'+mainC.chartID+'_'+i+'" style="color:'+mainC.datasets[i].params.color+'">'+mainC.datasets[i].params.name+'</FONT>';
				}
				if (numChkBoxes >= 5) {
					htmlString+='<LABEL><INPUT type="checkbox" id="datasetAutoExtend_'+mainC.chartID+'">Auto-extend</LABEL><BR>';
				}
				htmlString+=dsHtmlString;
				htmlString+="</FIELDSET>\n";
			}
			if (pointDatasets.length==0 && mainC.datasets.length==1 && (!mainC.datasets[0].params.type || (mainC.datasets[0].params.type === 'point' && !mainC.datasets[0].params.colorRule))) {pointDatasets.push(0);}
			
			for (let i=11; i<20; i++) { // custom elements [11..19]
				if (customFormElements[i]) {htmlString+=customFormElements[i];}
			}
			/* Selection buttons [index=20] */
			if (mainC.selectable || (customFormElements[20] && customFormElements[20][0]===true)) {
				let label=(customFormElements[20] && customFormElements[10][1])? customFormElements[20][1] : 'Selection';
				htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>'+label+':</B></LEGEND>';
				if (mainC.datasets.length > 1) {htmlString+='<INPUT type="checkbox" value="1" onclick="idvis.registeredCharts['+mainC.chartID+'].selectAllDatasets=this.checked"><FONT style="font-weight:bold;font-size:12px">Extend to all datasets</FONT><BR>';}
				htmlString+='<INPUT type="button" value="Select" style="font-weight:bold;font-size:12px;width:90px" onclick="idvis.lib.manageSelection(idvis.registeredCharts['+mainC.chartID+'],\'on\')">';
				htmlString+='<INPUT type="button" value="Un-Select" style="font-weight:bold;font-size:12px;width:90px" onclick="idvis.lib.manageSelection(idvis.registeredCharts['+mainC.chartID+'],\'off\')">';
				if (mainC.pointOnList) {
					if (typeof(mainC.pointOnList)=='function') { // default
						htmlString+='<BR><INPUT type="button" value="List selected" style="font-weight:bold;font-size:12px;width:180px" onclick="idvis.lib.listSelectedPoints(idvis.registeredCharts['+mainC.chartID+'],-1)">';
					}
					else { // array of ['action name',function] arrays
						for (let i=0; i<mainC.pointOnList.length; i++) {
							htmlString+='<BR><INPUT type="button" value="'+mainC.pointOnList[i][0]+' selected" style="font-weight:bold;font-size:12px;width:180px" onclick="idvis.lib.listSelectedPoints(idvis.registeredCharts['+mainC.chartID+'],'+i+')">';
						}
					}
				}
				htmlString+="</FIELDSET>\n";
			}
		}
		
		for (let i=21; i<30; i++) { // custom elements
			if (customFormElements[i]) {htmlString+=customFormElements[i];}
		}
		/* Search [index=30] */
		if (mainC.searchable || (customFormElements[30] && customFormElements[30][0]===true)) {
			let label=(customFormElements[30] && customFormElements[30][1])? customFormElements[30][1] : 'Search';
			htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>'+label+':</B></LEGEND>';
			let searchBoxID=mainC.divID+'_search',
				searchResDivID=mainC.divID+'_srchRes';
			htmlString+='<INPUT type="text" id="'+searchBoxID+'" style="width:140px" value=""><INPUT type="button" value="Go" style="font-weight:bold;font-size:12px;width:40px" onclick="idvis.lib.searchDataPoints(idvis.registeredCharts['+mainC.chartID+'],document.getElementById(\''+searchBoxID+'\').value,\''+searchResDivID+'\',';
			if (typeof(mainC.searchable)=='object') {
				let extSearchText=mainC.searchable.text || 'Extended search',
				extSearchChkID=mainC.divID+'_extSearch';
				htmlString+='\''+extSearchChkID+'\')"><BR><INPUT type="checkbox" id="'+extSearchChkID+'" value="1">'+extSearchText;
			}
			else {
				htmlString+='null)">';
			}
			htmlString+='<DIV id="'+searchResDivID+'"></DIV>';
			htmlString+="</FIELDSET>\n";
		}

		for (let i=31; i<40; i++) { // custom elements
			if (customFormElements[i]) {htmlString+=customFormElements[i];}
		}
		/* Editable threshold lines control [index=40] */
		if ((mainC.editableThresholds && mainC.editableThresholds.length) || mainC.dimArea || (customFormElements[40] && customFormElements[40][0]==true)) {
			let label=(customFormElements[40] && customFormElements[40][1])? customFormElements[40][1] : 'Thresholds';
			htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>'+label+':</B></LEGEND>';
			let editThres=false;
			if ((mainC.editableThresholds && mainC.editableThresholds.length) || mainC.dimArea || (customFormElements[40] && customFormElements[40][0]==true)) {
				let thresSelID=mainC.divID+'_thSel',
					thresBoxID=mainC.divID+'_thBox',
					dispValue='';
				if (mainC.editableThresholds.length > 1) {
					htmlString+='<SELECT id="'+thresSelID+'" style="font-size:12px;width:180px" onchange="idvis.lib.selectThreshold(idvis.registeredCharts['+mainC.chartID+'],this.value,\''+thresBoxID+'\')"><OPTION value="-1">-= Select =-</OPTION>';
					for (let i=0; i<mainC.editableThresholds.length; i++) {
						htmlString+='<OPTION value="'+i+'">'+mainC.editableThresholds[i].name+'</OPTION>';
					}
					htmlString+='</SELECT>';
				}
				else {
					htmlString+='<INPUT type="hidden" id="'+thresSelID+'" value="0"/><FONT style="font-size:12px;">'+mainC.editableThresholds[0].name+'</FONT>';
					dispValue=mainC.editableThresholds[0].initValue;
				}
				htmlString+='<INPUT type="text" id="'+thresBoxID+'" style="font-size:12px;width:60px" value="'+dispValue+'"><INPUT type="button" value="Apply" style="font-weight:bold;font-size:12px;width:60px" onclick="idvis.lib.updateThreshold(idvis.registeredCharts['+mainC.chartID+'],\''+thresSelID+'\',\''+thresBoxID+'\')">';
				editThres=true;
			}
			if (mainC.dimArea) {
				if (editThres) htmlString+='<BR>';
				htmlString+='<LABEL><INPUT type="checkbox" onclick="idvis.lib.activateDimArea(idvis.registeredCharts['+mainC.chartID+'],this.checked)"';
				if (mainC.dimArea.active===true) htmlString+=' checked';
				let label=mainC.dimArea.label || 'Dim Dim points based on thresholds';
				htmlString+='/>'+label+'</LABEL>';
			}
			htmlString+="</FIELDSET>\n";
		}

		for (let i=41; i<50; i++) { // custom elements
			if (customFormElements[i]) {htmlString+=customFormElements[i];}
		}
		/* Point highlighting [index=50] */
		if (mainC.allowHighlight || (customFormElements[50] && customFormElements[60][0]===true)) { // mainC.datasets.length &&
			let highlightDivID=mainC.divID+'_hlight';
			let label=(customFormElements[50] && customFormElements[50][1])? customFormElements[50][1] : 'Highlighting';
			htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>'+label+':</B></LEGEND><DIV id="'+mainC.divID+'_hlight">None</DIV></FIELDSET>\n';
		}

		for (let i=51; i<60; i++) { // custom elements
			if (customFormElements[i]) {htmlString+=customFormElements[i];}
		}

		/* Color management */
		if (!mainC.noColorEdition && (pointDatasets.length || (mainC.features && mainC.features.length))) {
			htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Color selection:</B></LEGEND>\n';
			htmlString+='<SELECT id="colorElementSel_'+mainC.chartID+'" onchange="idvis.lib.toggleColorSelection(idvis.registeredCharts['+mainC.chartID+'],this.value)"><OPTION value="none">-= Select =-</OPTION>';
			
			//Datasets
			if (pointDatasets.length) {
				htmlString+='<OPTGROUP label="'+datasetsLabel+':">';
				for (let i=0; i<pointDatasets.length; i++) {
					let idx=pointDatasets[i];
					htmlString+='<OPTION value="dataset:'+idx+'">'+mainC.datasets[idx].params.name+'</OPTION>';
				}
				htmlString+='</OPTGROUP>\n';
			}
			//Features
			if (mainC.features && mainC.features.length) {
				htmlString+='<OPTGROUP label="Features:">';
				for (let i=0; i<mainC.features.length; i++) {
					htmlString+='<OPTION value="feature:'+i+'">'+(mainC.features[i].name || idvis.ucFirst(mainC.features[i].type) + ' #'+i)+'</OPTION>';
				}
				htmlString+='</OPTGROUP>\n';
			}
			htmlString+='</SELECT>\n';
			//Color picker
			htmlString+='<DIV id="colorDIV_'+mainC.divID+'" style="display:none"><INPUT type="color" id="colorPicker_'+mainC.chartID+'" value="#000" style="height:20px;vertical-align:middle;">'; //  onchange="idvis.lib.changeElementColor(idvis.registeredCharts['+mainC.chartID+'],document.'+formName+'.colorElementSel.value,this.value)"
			htmlString+='<INPUT type="button" value="Apply" style="font-weight:bold;font-size:12px;width:60px;vertical-align:middle;" onclick="idvis.lib.changeElementColor(idvis.registeredCharts['+mainC.chartID+'],document.getElementById(\'colorElementSel_'+mainC.chartID+'\').value,document.getElementById(\'colorPicker_'+mainC.chartID+'\').value)"/></DIV>';
			htmlString+='</FIELDSET>\n';
		}

		/* Export image button [index=60] */
		if (mainC.exportAsImage || (customFormElements[60] && customFormElements[60][0]===true)) {
			if (htmlString) {htmlString+='<BR>';}
			if (mainC.exportAsImage[3]) { // image format
				htmlString+='<INPUT type="button" value="'+mainC.exportAsImage[0]+'" style="font-weight:bold;font-size:10px" onclick="idvis.lib.exportSVGtoImg(\''+mainC.getDivID()+'\',\''+mainC.exportAsImage[1]+'\',\''+mainC.exportAsImage[2]+'\',\''+mainC.exportAsImage[3]+'\')"/>\n';
			}
			else {
				htmlString+='<FONT style="font-weight:bold;font-size:11px">'+mainC.exportAsImage[0]+':</FONT>';
				//PNG
				htmlString+='<INPUT type="button" value="PNG" style="font-weight:bold;font-size:10px" onclick="idvis.lib.exportSVGtoImg(\''+mainC.getDivID()+'\',\''+mainC.exportAsImage[1]+'\',\''+mainC.exportAsImage[2]+'\',\'png\')"/>';
				//SVG
				htmlString+='<INPUT type="button" value="SVG" style="font-weight:bold;font-size:10px" onclick="idvis.lib.exportSVGtoImg(\''+mainC.getDivID()+'\',\''+mainC.exportAsImage[1]+'\',\''+mainC.exportAsImage[2]+'\',\'svg\')"/>\n';
			}
		}

		// htmlString+='</FORM>';
		mainC.formDiv.innerHTML=htmlString;

		/*Draw SVG icons for datasets if any */
		for (let i=0; i<datasetIcons.length; i++) {
			let [dsIdx,spanID]=datasetIcons[i];
			let c=Raphael(spanID,23,18);
//c.rect(0,0,23,18);
			if (!mainC.datasets[dsIdx].params.type || mainC.datasets[dsIdx].params.type === 'point') {
				mainC.datasets[dsIdx].params.icon=drawChartPoint(c,11,11,6,mainC.datasets[dsIdx].params.pointShape).attr({stroke:'none',fill:mainC.datasets[dsIdx].params.color});
			}
			else { // category box|bar|hist
				mainC.drawDatasetIcon(c,mainC.datasets[dsIdx]);
			}
		}

	};

	/******************* Chart drawing *********************/
	var initializeChart = function (mainC) { // *public*
		let paper=mainC.canvas,
			chartList=[];
		if (mainC.subChart) {
			for (let c=0; c<mainC.activeCharts.length; c++) {
				chartList.push(mainC.subChart[mainC.activeCharts[c]]);
			}
		}
		else {chartList.push(mainC);}
		
		var _handleDragContextDblClick=function() {
			return function() {if (this.data('chart').zoomable) {zoomIn(this.data('chart'));} else {_clearDragArea(this);}};
		};

		// Looping through displayed charts
		for (let c=0; c<chartList.length; c++) {
			let C=chartList[c],
//C.plotArea.attr({stroke:'#F00'});
				chartX=C.plotArea.attr('x')-0.5, // var chartX0=chartX-1;
				chartY=C.plotArea.attr('y')-1.5,
				chartW=C.plotArea.attr('width')+2,
				chartH=C.plotArea.attr('height')+2;

			/* Axes */
			C.usedAxes={};
			let axes=['X','Y','X2','Y2'];
			for (let i=0; i<axes.length; i++) {if (C['axis'+axes[i]+'text'] !== undefined) {C.usedAxes[axes[i]]=true;}}// axis is declared
//console.log(C.usedAxes);
			if (mainC.axisClosure) {paper.path('M'+chartX+' '+chartY+' l0 '+chartH+' l'+chartW+' 0 l0 -'+chartH+' Z');}
			else {
				if (C.usedAxes.X && !C.noAxisX) {paper.path('M'+chartX+' '+(chartY+chartH)+' l'+chartW+' 0');} // X
				if (C.usedAxes.Y && !C.noAxisY) {paper.path('M'+chartX+' '+chartY+' l0 '+chartH);} // Y
				if (C.usedAxes.X2 && !C.noAxisX2) {paper.path('M'+chartX+' '+chartY+' l'+chartW+' 0');} // X2
				if (C.usedAxes.Y2 && !C.noAxisY2) {paper.path('M'+(chartX+chartW)+' '+chartY+' l0 '+chartH);} // Y2
			}
			if (C.axisXtext) {
				let titleX=paper.text(chartX+chartW/2,chartY+chartH+30,C.axisXtext).attr({'font-weight':'bold','font-size':14});
				if (C.axisXtitle) C.axisXtitle=titleX; // used by pcaPlot.js to redraw title
			}
			if (C.axisYtext) {
				let tx=chartX-45, //Math.max(10,chartX-45);
					ty=chartY+chartH/2,
					titleY=paper.text(tx,ty,C.axisYtext).attr({'font-weight':'bold','font-size':14}).rotate(-90,tx,ty);
				if (C.axisYtitle) C.axisYtitle=titleY; // used by pcaPlot.js to redraw title
			}
			if (C.axisX2text) {
				let titleX2=paper.text(chartX+chartW/2,chartY-30,C.axisX2text).attr({'font-weight':'bold','font-size':14}); // Math.max(10,chartY-30)
				if (C.axisX2title) C.axisX2title=titleX2;
			}
			if (C.axisY2text) {
				let tx=chartX+chartW+45,
					ty=chartY+chartH/2,
					titleY2=paper.text(tx,ty,C.axisY2text).attr({'font-weight':'bold','font-size':14}).rotate(-90,tx,ty);
				if (C.axisY2title) C.axisY2title=titleY2;
			}


			/***** Drag selection area *****/
			if (C.dragContext) {
				C.dragArea=paper.rect(0,0,0,0,0)
				.attr({stroke:'#F00',fill:'#600','fill-opacity':0.1})
				.data({startX:0,startY:0,status:'off',chart:C})
				//.dblclick(function() {if (this.data('chart').zoomable) {zoomIn(this.data('chart'));} else {_clearDragArea(this);}})
				.dblclick(_handleDragContextDblClick())
				.hide();

				if (C.zoomable) {
					if (!C.hideZoom) {
						/* Zoom info */
						C.zoomText=paper.text(chartX,10,'').attr({'font-weight':'bold','font-size':12,'text-anchor':'start'});
					}
					/* Drag area dblclick */
					//C.dragArea=paper.rect(0,0,0,0,0)
					//.attr({stroke:'#F00',fill:'#600','fill-opacity':0.1})
					//.data({startX:0,startY:0,status:'off',chart:C})
					//C.dragArea.dblclick(function() {zoomIn(this.data('chart'));});
					//.hide();
				}
			}
		}

		/***** Data points to SVG points *****/
		if (mainC.chartType !== 'categoryPlot' && mainC.datasets) {
			var _addHoverToPoint=function(p) {
				p.hover(function(e){setLabelDisplay(this,'on','max',e);},function(){setLabelDisplay(this,'off');})
			}
			//var _handlePointHover=function(type) {
			//	return (type==='on')? function() {setLabelDisplay(this,'on','max');} : function() {setLabelDisplay(this,'off');};
			//};
			var _handlePointClick=function() {
				return function() {if (idvis.modKeyPressed && this.data('js').dataset.params.chart.onPointExclusion) {setPointExclusion(this);} else {setPointSelection(this,'auto','max');}};
			};
			for (let dsIdx=0; dsIdx < mainC.datasets.length; dsIdx++) {
				const pointSizeRule=(mainC.datasets[dsIdx].params.pointSizeRule)? mainC.datasets[dsIdx].params.pointSizeRule : {min:1,max:20,scale:2};
				const pointShape=mainC.datasets[dsIdx].params.pointShape;
				for (let i=0; i < mainC.datasets[dsIdx].data.length; i++) {
					let dp=mainC.datasets[dsIdx].data[i],
					//var usedSize=(dp.size)? Math.min(dp.size,20) : 5;
						usedSize=(dp.size)? Math.max(Math.min(dp.size,pointSizeRule.max),pointSizeRule.min) : 5, // max point size
					//var p=paper.circle(-50,-50,2*Math.sqrt(usedSize/3.14))
					//p=paper.circle(-50,-50,pointSizeRule.scale*Math.sqrt(usedSize/3.14))
					p=drawChartPoint(paper,-50,-50,pointSizeRule.scale*Math.sqrt(usedSize/3.14),pointShape)
						.attr({stroke:'none',fill:mainC.datasets[dsIdx].params.color,'fill-opacity':mainC.pointOpacity}) //,'fill-opacity':0.65   title:'point '+i,cursor:"crosshair"
						.data({js:dp,showLabel:null}) //x:pX,y:pY,
						//.hover(function(){setLabelDisplay(this,'on','max');},function(){setLabelDisplay(this,'off');})
						//.click(function(){if (idvis.modKeyPressed && this.data('js').dataset.params.chart.onPointExclusion) {setPointExclusion(this);} else {setPointSelection(this,'auto','max');}});
						//.hover(_handlePointHover('on'),_handlePointHover('off'))
						.click(_handlePointClick());
						_addHoverToPoint(p);
					dp.point=p;
					if (dp.valueList) { // multi-point data
						dp.pointList=[];
						//var r=Math.round(1.4*Math.sqrt(usedSize/3.14)); // smaller than mater point
						let r=(mainC.datasets[dsIdx].params.multiValuePoint.size)? 2*Math.sqrt(mainC.datasets[dsIdx].params.multiValuePoint.size/3.14) : 1.4*Math.sqrt(usedSize/3.14); // smaller than mater point
						if (mainC.datasets[dsIdx].params.multiValuePoint.linked) {
							for (let i=0; i<dp.valueList.x.length; i++) { // only 1 point per x,y pair
								// dp.pointList.push(paper.circle(-50,-50,r).attr({stroke:'none',fill:mainC.datasets[dsIdx].params.color,'fill-opacity':mainC.pointOpacity}).hide()); // hidden
								dp.pointList.push(drawChartPoint(paper,-50,-50,r,pointShape).attr({stroke:'none',fill:mainC.datasets[dsIdx].params.color,'fill-opacity':mainC.pointOpacity}).hide()); // hidden
							}
						}
						else {
							for (let axis in dp.valueList) {
								if (dp.valueList.hasOwnProperty(axis)) {
									for (let j=0; j<dp.valueList[axis].length; j++) {
										// dp.pointList.push(paper.circle(-50,-50,r).attr({stroke:'none',fill:mainC.datasets[dsIdx].params.color,'fill-opacity':mainC.pointOpacity}).hide()); // hidden
										dp.pointList.push(drawChartPoint(paper,-50,-50,r,pointShape).attr({stroke:'none',fill:mainC.datasets[dsIdx].params.color,'fill-opacity':mainC.pointOpacity}).hide()); // hidden
									}
								}
							}
						}
						p.toFront();
					}
				}
			}
		}

		/***** Highlighting legends *****/
		mainC.highlightLegends=paper.set();

		/***** Chart-specific initialization *****/
		if (mainC.initializeChart) {
			mainC.initializeChart();
		}

		// Looping through displayed charts
		for (let c=0; c<chartList.length; c++) {
			let C=chartList[c];
			/***** Compute default range *****/
			computeDefaultRange(C);

			/***** Plot chart *****/
			plotChart(C);
		}
//,minRangeX,maxRangeX,optRangeX[0],optRangeX[2],minRangeY,maxRangeY,optRangeY[0],optRangeY[2]);
	};

	var computeDefaultRange = function(C) { // *public*
		let chartW=C.plotArea.attr('width'),
			chartH=C.plotArea.attr('height'),
			ref=C.chartSettings.reference,
			minValues={},
			maxValues={};
		if (C.sameScale) { // same scale for both axes at start (PCA)
			let maxRange=0,
				ranges={};
			for (let axis in C.usedAxes) {
				if (C.usedAxes.hasOwnProperty(axis)) {
					ranges[axis]=C['maxValue'+axis]-C['minValue'+axis];
					if (maxRange < ranges[axis]) maxRange=ranges[axis];
				}
			}
			for (let axis in ranges) {
				if (ranges.hasOwnProperty(axis)) {
					let extraRange=(maxRange-ranges[axis])/2;
					minValues[axis]=C['minValue'+axis]-extraRange;
					maxValues[axis]=C['maxValue'+axis]+extraRange;
				}
			}
		}
		else {
			for (let axis in C.usedAxes) {
				if (C.usedAxes.hasOwnProperty(axis)) {
					minValues[axis]=C['minValue'+axis];
					maxValues[axis]=C['maxValue'+axis];
				}
			}
		}

		for (let axis in C.usedAxes) {
			if (C.usedAxes.hasOwnProperty(axis)) {
				let axisSize=(axis.match('X'))? chartW : chartH,
					optRange=getChartScaleRange(true,minValues[axis],maxValues[axis],axisSize);
				ref[axis]={};
				ref[axis].minRange=(!C['force'+axis+'to0'])? optRange[0] : (C['force'+axis+'to0']==1)? 0 : Math.max(optRange[0],-optRange[3]/5); // 0 or null: auto, 1: 0, 2: ~0
				ref[axis].optMinRange=optRange[1];
				ref[axis].maxRange=optRange[2];
				ref[axis].tickSize=optRange[3];
				ref[axis].pix2valRatio=(ref[axis].maxRange-ref[axis].minRange)/axisSize;
			}
		}

		// set current=reference
		for (let axis in ref) {
			if (ref.hasOwnProperty(axis)) {
				C.chartSettings.current[axis]={};
				for (let prop in ref[axis]) {
					if (ref[axis].hasOwnProperty(prop)) {
						C.chartSettings.current[axis][prop]=ref[axis][prop];
					}
				}
			}
		}
		ref.zoomX=ref.zoomY=C.chartSettings.current.zoomX=C.chartSettings.current.zoomY=1;
	};

	var plotChart = function(C) { // *public*
		let mainC=(C.mainChart)? C.mainChart : C,
			paper=mainC.canvas,
			chartX=C.plotArea.attr('x'), chartX0=chartX-0.5,
			chartY=C.plotArea.attr('y'), chartY0=chartY-1.5,
			chartW=C.plotArea.attr('width'),
			chartH=C.plotArea.attr('height'),
			settings=C.chartSettings.current;

		/***** Thresold lines *****/
		//if (C.thresholdLines) {
		//	for (let i=0; i<C.thresholdLines.length; i++) {
		//		_moveThresholdLine(C.thresholdLines[i]);
		//	}
		//}


		/***** Features *****/
		if (C.features) {
			for (let i=0; i<C.features.length; i++) {
				let f=C.features[i],
					anchorPoints=[];
//console.log(i+': '+f.type);
				/* Threshold line */
				if (f.type==='line') {
					if (f.axis.match('X')) {
						f.pathPos=chartX0+Math.round((f.value-settings[f.axis].minRange)/settings[f.axis].pix2valRatio);
						if (f.pathPos >= chartX && f.pathPos <= chartX0+chartW) {
							anchorPoints=[[f.pathPos,chartY0],[f.pathPos,chartY+chartH]];
						}
					}
					else {
						f.pathPos=chartY0+chartH-Math.round((f.value-settings[f.axis].minRange)/settings[f.axis].pix2valRatio);
						if (f.pathPos >= chartY && f.pathPos <= chartY0+chartH) {
							anchorPoints=[[chartX0,f.pathPos],[chartX+chartW,f.pathPos]];
						}
						//f.pathPos=f.pathPos;
					}
				}
				/* Function */
				else if (f.type==='function') {
//console.log(f);
					//let axisY=(f.axis==='X')? 'Y' : 'Y2',
					const axisXdelta=(settings[f.axis].maxRange-settings[f.axis].minRange)/20;
					for (let x0=settings[f.axis].minRange-2*axisXdelta; x0 <= settings[f.axis].maxRange+2*axisXdelta; x0+=axisXdelta) {
						let x=chartX0+Math.round((x0-settings[f.axis].minRange)/settings[f.axis].pix2valRatio),
							y=chartY0+chartH-Math.round((f.data(x0)-settings[f.axisY].minRange)/settings[f.axisY].pix2valRatio);
						anchorPoints.push([x,y]);
					}
//console.log(anchorPoints);
				}
				/* Custom */
				else if (f.type==='custom') {
					//let axisY=(f.axis==='X')? 'Y' : 'Y2';
					for (let d=0; d<f.data.length; d++) {
						let x=chartX0+Math.round((f.data[d][0]-settings[f.axis].minRange)/settings[f.axis].pix2valRatio),
							y=chartY0+chartH-Math.round((f.data[d][1]-settings[f.axisY].minRange)/settings[f.axisY].pix2valRatio);
						anchorPoints.push([x,y]); // TODO: filter for visible section after zoom
					}
				}
//console.log(anchorPoints);
				/* Text */
				else if (f.type==='text') {
					//let axisY=(f.axis==='X')? 'Y' : 'Y2',
					let x=chartX0+Math.round((f.x-settings[f.axis].minRange)/settings[f.axis].pix2valRatio),
						y=chartY0+chartH-Math.round((f.y-settings[f.axisY].minRange)/settings[f.axisY].pix2valRatio);
					f.path=paper.text(x,y,f.label).attr({'text-anchor':f.anchor,fill:f.color,'font-size':f.size,'font-weight':f.weight});
				}

				drawFeature(C,f,anchorPoints);
/*				
				if (anchorPoints.length >= 2) {
					f.path=drawPath(C,f,anchorPoints);
					f.path.data({js:f});
					if (f.label) {
						f.path.hover(
							function(e) {
								this.attr({'stroke-width':3});
								let ol=this.data('js');
								//ol.popup=drawLabel(paper,this,ol.chart.plotArea.attr('x'),ol.chart.plotArea.attr('y')+ol.chart.plotArea.attr('height'),1,1,ol.label,ol.color);
								ol.popup=drawLabel(paper,this,e.layerX,e.layerY,1,1,ol.label,ol.color);
							},
							function(){
								this.attr({'stroke-width':1});
								let ol=this.data('js');
								ol.popup.remove();
								ol.popup=null;
							}
						);
					}
				}
*/
			}
		}


		/***** Axes ticks *****/
		let posX,posY;
		for (let axis in C.usedAxes) {
			if (!C.usedAxes.hasOwnProperty(axis)) {continue;}
			// X,X2 ticks
			if (axis.match('X')) {
				if (C['noTicks'+axis]) {
					if (C['noTicks'+axis] !== true) {
						posX=chartX0+Math.round(chartW/2);
						posY=(axis=='X')? chartY0+chartH+10 : chartY0-10;
						paper.text(posX,posY,C.noTicksX).attr({'font-weight':'bold','font-size':12});
					}
				}
				else {
					let sign;
					if (axis=='X') {sign=1; posY=chartY+chartH+1;}
					else {sign=-1; posY=chartY-1;} // posY 2 px below/above chart
					let tickSizeX=settings[axis].tickSize,
						fixedX=(tickSizeX < 1)? (Math.round(1/tickSizeX)+'').length : (Math.round(tickSizeX)==tickSizeX)? 0 : 1,
						tickVx=settings[axis].optMinRange;
					while (tickVx <= settings[axis].maxRange) {
						posX=chartX0+Math.round((tickVx-settings[axis].minRange)/settings[axis].pix2valRatio);
						if (posX >= chartX0) {
							C.chartMarks.push(paper.path('M'+posX+' '+posY+' l0 '+(sign*5)));
							C.chartMarks.push(paper.text(posX,posY+(sign*10),tickVx.toFixed(fixedX)));
						}
						tickVx+=tickSizeX;
					}
				}
			}
			// Y,Y2 ticks
			if (axis.match('Y')) {
				if (C['noTicks'+axis]) {
					if (C['noTicks'+axis] !== true) {
						posX=(axis=='Y')? chartX0-10 : chartX0+chartW+10;
						posY=chartY0+Math.round(chartH/2);
						paper.text(posX,posY,C.noTicksY).attr({'font-weight':'bold','font-size':12}).rotate(-90,posX,posY);
					}
				}
				else {
					let sign,txtAnch;
					if (axis=='Y') {sign=-1; posX=chartX0; txtAnch='end';}
					else {sign=1; posX=chartX+chartW+1; txtAnch='start';}
					let tickSizeY=settings[axis].tickSize,
						fixedY=(tickSizeY < 1)? (Math.round(1/tickSizeY)+'').length : (Math.round(tickSizeY)==tickSizeY)? 0 : 1,
						tickVy=settings[axis].optMinRange;
					while (tickVy <= settings[axis].maxRange) {
						posY=chartY0+chartH-Math.round((tickVy-settings[axis].minRange)/settings[axis].pix2valRatio);
						if (posY <= chartY0+chartH) {
							C.chartMarks.push(paper.path('M'+posX+' '+posY+' l'+(sign*5)+' 0'));
							C.chartMarks.push(paper.text(posX+(sign*6),posY,tickVy.toFixed(fixedY)).attr({'text-anchor':txtAnch}));
						}
						tickVy+=tickSizeY;
					}
				}
			}
		}
//console.log('PresX='+precisX+' ('+tickSizeX+'), PresY='+precisY+' ('+tickSizeY+')');

		//Zooms
		if (C.zoomable && !C.hideZoom) {C.zoomText.attr('text','ZoomX: '+settings.zoomX.toPrecision(2)+' - ZoomY: '+settings.zoomY.toPrecision(2));}

		/***** Ploting data points *****/
		let existLines=false;
		if (mainC.datasets) {
			var _handleLineHover=function(type) {
				if (type==='on') {
					return function(e) {
						this.attr({'stroke-width':Math.max(this.attr('stroke-width')*2,idvis.pathHoverWidth)});
						let ol=this.data('js'),
							paper=(C.mainChart)? C.mainChart.canvas : C.canvas;
						//ol.popup=drawLabel(paper,this,ol.chart.plotArea.attr('x'),ol.chart.plotArea.attr('y')+ol.chart.plotArea.attr('height'),1,1,ol.label,ol.color);
						ol.popup=drawLabel(paper,this,e.layerX,e.layerY,1,1,ol.label,ol.color);
					};
				}
				else {
					return function(){
						let ol=this.data('js');
						this.attr({'stroke-width':ol.width});
						ol.popup.remove();
						ol.popup=null;
					};
				}
			};
			let datasets=mainC.datasets;
			for (let dsIdx=0; dsIdx < datasets.length; dsIdx++) {
				//var pathStrg='',startLine=true;
				let dsX=(datasets[dsIdx].params.axisX)? datasets[dsIdx].params.axisX : 'X',
					dsY=(datasets[dsIdx].params.axisY)? datasets[dsIdx].params.axisY : 'Y',
					pathPoints=[];
				//var pointOut=false; var pointIn=false;
				for (let i=0; i < datasets[dsIdx].data.length; i++) {
					let dp=datasets[dsIdx].data[i];
					if (dp.subChart && dp.subChart !== C) continue;
					const pointShape=mainC.datasets[dsIdx].params.pointShape;
					let p=dp.point;
					// point in range
					if (dp.getX() >= settings[dsX].minRange && dp.getX() <= settings[dsX].maxRange && dp.getY() >= settings[dsY].minRange && dp.getY() <= settings[dsY].maxRange) {
						pointIn=true;
						let pX=Math.round((dp.getX()-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0,
							pY=chartY0+chartH-Math.round((dp.getY()-settings[dsY].minRange)/settings[dsY].pix2valRatio);
						// p.attr({cx:pX,cy:pY});
						moveChartPointTo(p,pX,pY,pointShape);
						if (datasets[dsIdx].params.visible && !dp.isHidden) {
							p.show();
							if ((mainC.connectPoints || datasets[dsIdx].dataConnection) && !dp.noLine) {
								/*
								//if (pointOut) {
								//	var prevDp=datasets[dsIdx].data[i-1];
								//	var prevX=Math.round((prevDp.getX()-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0;
								//	var prevY=chartY0+chartH-Math.round((prevDp.getY()-settings[dsY].minRange)/settings[dsY].pix2valRatio);
								//	prevDp.point.attr({cx:prevX,cy:prevY});
								//	pathPoints.push(prevDp.point);
								//	pointOut=false;
								//}
								*/
								pathPoints.push(p);
							}
							if (p.data('showLabel')) { // showLabel can be on but not labelObj due to panning
								if (C.dragContext=='pan' && p.data('labelObj')) {
									let tSet=p.data('labelObj');
									tSet.translate(pX-tSet[0].data('x'),pY-tSet[0].data('y'));
									tSet[0].data({x:pX,y:pY});
									tSet.show();
								}
								else {
									_emphasizePoint(p,'on',p.data('showLabel'),pX,pY);
								}
							}
							if (dp.pointList) { // multi-point data
								if (datasets[dsIdx].params.multiValuePoint.linked) { // 1 point for each x,y pair
									for (let j=0; j<dp.valueList.x.length; j++) {
										let pX2=Math.round((dp.valueList.x[j]-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0,
											pY2=chartY0+chartH-Math.round((dp.valueList.y[j]-settings[dsY].minRange)/settings[dsY].pix2valRatio);
										// dp.pointList[j].attr({cx:pX2,cy:pY2});
										moveChartPointTo(dp.pointList[j],pX2,pY2,pointShape);
									}
								}
								else {
									let k=-1;
									for (let axis in dp.valueList) {
										if (dp.pointList.hasOwnProperty(axis)) {
											let pX2,pY2;
											if (axis=='x') {pY2=pY;} else {pX2=pX;}
											for (let j=0; j<dp.valueList[axis].length; j++) {
												if (axis=='x') {pX2=Math.round((dp.valueList.x[j]-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0;}
												else {pY2=chartY0+chartH-Math.round((dp.valueList.y[j]-settings[dsY].minRange)/settings[dsY].pix2valRatio);}
												k++;
												// dp.pointList[k].attr({cx:pX2,cy:pY2});
												moveChartPointTo(dp.pointList[k],pX2,pY2,pointShape);
											}
										}
									}
								}
							}
						}
					}
					// point out of range
					else {
						if (p.data('labelObj') && datasets[dsIdx].params.visible && C.dragContext=='pan') {
							p.data('labelObj').hide();
						}
						//pointOut=true;
						// Line => compute out of range point & add to path list
						if ((mainC.connectPoints || datasets[dsIdx].dataConnection) && !dp.noLine) { //pointIn &&
							let pX=Math.round((dp.getX()-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0,
								pY=chartY0+chartH-Math.round((dp.getY()-settings[dsY].minRange)/settings[dsY].pix2valRatio);
							// p.attr({cx:pX,cy:pY});
							moveChartPointTo(p,pX,pY,pointShape);
							pathPoints.push(p);
							//pointIn=false;
						}
					}
				}
				var isLinearRegression=false;
				if (datasets[dsIdx].dataConnection && datasets[dsIdx].dataConnection.lineType==='regression') {
					var lr;
					if (datasets[dsIdx].dataConnection.lr) {
						lr=datasets[dsIdx].dataConnection.lr;
					}
					else { // 1srt draw: compute regression
						//Build x,y arrays
						let xValues=[], yValues=[];
						for (let i=0; i<pathPoints.length; i++) {
							xValues.push(pathPoints[i].data('js').getX());
							yValues.push(pathPoints[i].data('js').getY());
						}
						//Compute linear regression parameters
						lr=linearRegression(xValues,yValues);
						datasets[dsIdx].dataConnection.lr=lr;
						datasets[dsIdx].dataConnection.label=datasets[dsIdx].params.name+'\nSlope: '+idvis.formatNumber(lr.slope)+'\nIntercept: '+idvis.formatNumber(lr.intercept)+'\nR2: '+idvis.formatNumber(lr.r2);
//console.log(lr,datasets[dsIdx].dataConnection.label); //slope, intercept, r, r2
					}
					let x0=chartX0,
						x1=chartX0+chartW,
						y0=chartY0+chartH-Math.round((((lr.slope * settings[dsX].minRange) + lr.intercept)-settings[dsY].minRange)/settings[dsY].pix2valRatio);
						y1=chartY0+chartH-Math.round((((lr.slope * settings[dsX].maxRange) + lr.intercept)-settings[dsY].minRange)/settings[dsY].pix2valRatio);
					pathPoints=[[x0,y0],[x1,y1]];  // overwrite pathPoints
					isLinearRegression=true;
//console.log(pathPoints);	
				}
						
				if (pathPoints.length) {
//console.log(dsIdx+': '+pathPoints.length);					
					datasets[dsIdx].dataConnection.svg=drawPath(C,datasets[dsIdx],pathPoints).toBack(); // move behind points
					existLines=true;
					
					if (isLinearRegression) {
						let dsLine=datasets[dsIdx].dataConnection;
						dsLine.svg.data({js:dsLine});
						dsLine.svg.hover(_handleLineHover('on'),_handleLineHover('off'));
						/*
						dsLine.svg.hover(
							function(e) {
								this.attr({'stroke-width':Math.max(this.attr('stroke-width')*2,idvis.pathHoverWidth)});
								let ol=this.data('js'),
									paper=(C.mainChart)? C.mainChart.canvas : C.canvas;
								//ol.popup=drawLabel(paper,this,ol.chart.plotArea.attr('x'),ol.chart.plotArea.attr('y')+ol.chart.plotArea.attr('height'),1,1,ol.label,ol.color);
								ol.popup=drawLabel(paper,this,e.layerX,e.layerY,1,1,ol.label,ol.color);
							},
							function(){
								let ol=this.data('js');
								this.attr({'stroke-width':ol.width});
								ol.popup.remove();
								ol.popup=null;
							}
						);
						*/
					}
				}
			}
		}
		if (existLines) { // move chart panels behind line paths
			C.plotArea.toBack();
			mainC.backPanel.toBack();
		}

		if (C.plotChart) { // Chart-specific plot
			C.plotChart();
		}
	};
	
	/***** Feature *****/
	var drawFeature = function(C,f,anchorPoints) {
		if (anchorPoints.length >= 2) {
			f.path=drawPath(C,f,anchorPoints);
			f.path.data({js:f});
			if (f.label) {
				f.path.hover(
					function(e) {
						this.attr({'stroke-width':Math.max(this.attr('stroke-width')*2,idvis.pathHoverWidth)});
						let ol=this.data('js'),
						    paper=(C.mainChart)? C.mainChart.canvas : C.canvas;
						//ol.popup=drawLabel(paper,this,ol.chart.plotArea.attr('x'),ol.chart.plotArea.attr('y')+ol.chart.plotArea.attr('height'),1,1,ol.label,ol.color);
						ol.popup=drawLabel(paper,this,e.layerX,e.layerY,1,1,ol.label,ol.color);
					},
					function(){
						let ol=this.data('js');
						this.attr({'stroke-width':ol.width});
						ol.popup.remove();
						ol.popup=null;
					}
				);
			}
		}
	};

	/*********** Point shape handling functions ************/
	const datasetIcons=['circle','square','triangle','star','diamond','triangle_down','plus','triangle_right','cross','triangle_left','star6'];
	var generateDatasetAES = function(mainC,dsIdx,set) {
		if (!mainC.multiPointDatasetRule) {mainC.multiPointDatasetRule='pointShape';} // What to change first when adding new dataset: pointShape, color
		if (!mainC.currentAESindexes) {
			mainC.currentAESindexes={pointShape:-1,color:-1}; // pointShape, color
			mainC.usedAES=[];
		}
		let dsPointShape,dsColor;
		if (set.pointShape && set.color) {
			[dsPointShape,dsColor]=[_decodePointShape(set.pointShape),set.color];
			let matchedDsIdx=mainC.usedAES.indexOf(dsPointShape+':'+dsColor);
			if (matchedDsIdx >= 0) { // a dataset with same AES exists => update it
//console.log('MATCH_USER',matchedDsIdx);
				[mainC.datasets[matchedDsIdx].params.pointShape,mainC.datasets[matchedDsIdx].params.color]=evaluateAES(mainC,dsPointShape,dsColor);
			}
		}
		else if (set.pointShape) { // compute dsColor
			let colorIdx=Math.max(0,mainC.currentAESindexes.color);
			[dsPointShape,dsColor]=evaluateAES(mainC,_decodePointShape(set.pointShape),idvis.colorList[colorIdx],'pointShape');
		}
		else if (set.color) { // compute dsPointShape
			let shapeIdx=Math.max(0,mainC.currentAESindexes.pointShape);
			[dsPointShape,dsColor]=evaluateAES(mainC,datasetIcons[shapeIdx],set.color,'color');
		}
		else { // compute both
			let [aes1,aes2,param1,param2]=(mainC.multiPointDatasetRule==='color')? [idvis.colorList,datasetIcons,'color','pointShape'] : [datasetIcons,idvis.colorList,'pointShape','color'];
			mainC.currentAESindexes[param1]++;
			if (mainC.currentAESindexes[param1] === 0) {mainC.currentAESindexes[param2]=0;} // initialize 2nd dim aes
			if (mainC.currentAESindexes[param1] === aes1.length) {
				mainC.currentAESindexes[param1]=0;
				mainC.currentAESindexes[param2]++;
				if (mainC.currentAESindexes[param2] === aes2.length) {
					mainC.currentAESindexes[param2]=0;
				}
			}
			[dsPointShape,dsColor]=evaluateAES(mainC,datasetIcons[mainC.currentAESindexes.pointShape],idvis.colorList[mainC.currentAESindexes.color]);
		}
		mainC.usedAES.push(dsPointShape+':'+dsColor);
//console.log(dsIdx,mainC.usedAES);
		return [dsPointShape,dsColor];
	};
	
	var evaluateAES=function(mainC,dsPointShape,dsColor,lockedAes='-') {
		let [aes1,aes2,param1,param2]=(mainC.multiPointDatasetRule==='color')? [idvis.colorList,datasetIcons,'color','pointShape'] : [datasetIcons,idvis.colorList,'pointShape','color'],
			[newPointShape,newColor]=[dsPointShape,dsColor], // default
			maxLoopCount=(lockedAes===param1)? aes2.length : (lockedAes===param2)? aes1.length : aes1.length*aes2.length,
			loopCount=0;
		while (mainC.usedAES.indexOf(newPointShape+':'+newColor) >= 0) { // a dataset with same AES exists => find new aes pair
//console.log('MATCH_LOOP',lockedAes);
			loopCount++;
			if (lockedAes===param1) {
				mainC.currentAESindexes[param2]++;
				mainC.currentAESindexes[param2] %= aes2.length;
			}
			else if (lockedAes===param2) {
				mainC.currentAESindexes[param1]++;
				mainC.currentAESindexes[param1] %= aes1.length;
			}
			else { // no locked AES => loop on both
				mainC.currentAESindexes[param1]++;
				if (mainC.currentAESindexes[param1] >= aes1.length) { // Entire first aes tried, change 2nd
					mainC.currentAESindexes[param1]=0;
					mainC.currentAESindexes[param2]++;
					mainC.currentAESindexes[param2] %= aes2.length;
				}
			}
			newPointShape=datasetIcons[mainC.currentAESindexes.pointShape];
			newColor=idvis.colorList[mainC.currentAESindexes.color];
			if (loopCount === maxLoopCount) {break;}
		}
		return [newPointShape,newColor];
	};

	function _decodePointShape(shapeString='circle') { // public default circle
		shapeString=shapeString.toLowerCase();
		let pointShape=
			(idvis.isNumber(shapeString))? datasetIcons[shapeString % datasetIcons.length] :
			(datasetIcons.indexOf(shapeString) >= 0)? shapeString : // full name found
			(shapeString.match('^st6'))? 'star6' :
			(shapeString.match('^st'))? 'star' :
			(shapeString.match('^s'))? 'square' :
			(shapeString.match('^d'))? 'diamond' :
			(shapeString.match('^td'))? 'triangle_down' :
			(shapeString.match('^tr'))? 'triangle_right' :
			(shapeString.match('^tl'))? 'triangle_left' :
			(shapeString.match('^t'))? 'triangle' :
			(shapeString==='+' || shapeString.match('^p'))? 'plus' :
			(shapeString==='x' || shapeString.match('^cr'))? 'cross' :
			'circle';
//console.log(shapeString,'->',pointShape);
		return pointShape;
	}
	var drawChartPoint=function(paper,x,y,radius,shape='circle') { // public
		//shape='star';
		let p;
		if (shape==='circle') {
			p=paper.circle(x,y,radius); // (S=PI * r2)
		}
		else if (shape.match('square|diamond')) {
			let midWidth=radius/1.27, // scale square vs circle PI/4 (S=4 * r2)
				width=2*midWidth;
			p=paper.rect(x-midWidth,y-midWidth,width,width);
			if (shape==='diamond') {
				p.rotate(45);
			}
		}
		else if (shape.match('triangle|star6')) { // triangle path
			let c2top=radius*1.3, // scale triangle vs circle radius: PI/1.3 (S=1.3 * r2)
				height=c2top*1.5, // vertical height of triangle (pointing up)
				c2base=c2top/2,
				midSide=height/Math.sqrt(3),
				side=2*midSide; // length of triangle side
				//minY=y-c2top,
				//maxY=y+c2base; //midWidth;	
			if (shape==='triangle_down') {
				//p=paper.path('M'+x+' '+y+'M'+x+' '+maxY+'l-'+midSide+' -'+height+'l'+side+' 0L'+x+' '+maxY);
				p=paper.path('M'+x+' '+y+'m0 '+c2top+'l-'+midSide+' -'+height+'l'+side+' 0l-'+midSide+' '+height);
			}
			else if (shape==='triangle_right') {
				p=paper.path('M'+x+' '+y+'m'+c2top+' 0l-'+height+' '+midSide+'l0 -'+side+'l'+height+' '+midSide);
			}
			else if (shape==='triangle_left') {
				p=paper.path('M'+x+' '+y+'m-'+c2top+' 0l'+height+' -'+midSide+'l0 '+side+'l-'+height+' -'+midSide);
			}
			else if (shape==='star6') {
				//p=paper.path('M'+x+' '+y+'M'+x+' '+minY+'l'+midSide+' '+height+'l-'+side+' 0L'+x+' '+minY+'M'+x+' '+maxY+'l-'+midSide+' -'+height+'l'+side+' 0M'+x+' '+maxY);
				//p=paper.path('M'+x+' '+y+'m0 -'+c2top+'l'+midSide+' '+height+'l-'+side+' 0l'+midSide+' -'+height+'M'+x+' '+y+'m0 '+c2top+'l-'+midSide+' -'+height+'l'+side+' 0l-'+midSide+' '+height);
				p=paper.path('M'+x+' '+y+'m0 -'+c2top+'l'+midSide+' '+height+'l-'+side+' 0l'+midSide+' -'+height+'m0 '+(height+c2base)+'l-'+midSide+' -'+height+'l'+side+' 0l-'+midSide+' '+height);
			}
			else { // triangle up
				//p=paper.path('M'+x+' '+y+'M'+x+' '+minY+'l'+midSide+' '+height+'l-'+side+' 0L'+x+' '+minY);
				p=paper.path('M'+x+' '+y+'m0 -'+c2top+'l'+midSide+' '+height+'l-'+side+' 0l'+midSide+' -'+height);
			}
		}
		else if (shape==='star') { // sin: X axis, cos: Y axis
			let width=radius*1.2,
				sinW72=width*0.951, // sin(72°)
				cosW72=width*0.31, // cos(72°)
				sinW144=width*0.588, // sin(144°)
				cosW144=width*0.81; // -cos(144°)
			p=paper.path('M'+x+' '+y+'m0 -'+width+'l'+sinW144+' '+(width+cosW144)+'l-'+(sinW144+sinW72)+' -'+(cosW144+cosW72)+'l'+(sinW72+sinW72)+' 0l-'+(sinW72+sinW144)+' '+(cosW72+cosW144)+'l'+sinW144+' -'+(cosW144+width));
		}
		else if (shape==='plus') {
			let height=radius*2,
				//thick=Math.max(1,Math.round(radius/2.5)),
				thick=Math.max(1,radius/2.5),
				midThick=thick/2;
			p=paper.path('M'+x+' '+y+'m-'+midThick+' -'+radius+'l'+thick+' 0l0 '+height+'l-'+thick+' 0l0 -'+height+'M'+(x-radius)+' '+(y-midThick)+'l'+height+' 0l0 '+thick+'l-'+height+' 0l0 -'+thick);
		}
		else if (shape==='cross') {
			let dMidWidth=radius*0.71, // 0.71 <- sin/cos(45°)
				dWidth=dMidWidth*2,
				dThick=Math.max(1,dMidWidth/2.5),
				dMidThick=dThick/2,
				lMidWidth=dMidWidth-dMidThick,
				hMidWidth=dMidWidth+dMidThick;
			p=paper.path('M'+x+' '+y+'m'+lMidWidth+' -'+hMidWidth+'l'+dThick+' '+dThick+'l-'+dWidth+' '+dWidth+'l-'+dThick+' -'+dThick+'l'+dWidth+' -'+dWidth
			+'M'+(x-hMidWidth)+' '+(y-lMidWidth)+'l'+dThick+' -'+dThick+'l'+dWidth+' '+dWidth+'l-'+dThick+' '+dThick+'l-'+dWidth+' -'+dWidth
			);
		}
		return p;
	};
	var moveChartPointTo=function(p,x,y,shape) {
		if (shape.match('circle|square|diamond')) {
			let midWidth;
			switch (shape) {
				case 'square':
					midWidth=p.attr('width')/2;
					p.attr({x:x-midWidth,y:y-midWidth});
					break;
				case 'diamond':
					midWidth=p.attr('width')/2;
					p.rotate(-45).attr({x:x-midWidth,y:y-midWidth}).rotate(45);
					break;
				default: // circle
					p.attr({cx:x,cy:y});
			}
		}
		else { // path
			x-=p.attr('path')[0][1]; // Raphael.transform() is relative
			y-=p.attr('path')[0][2];
			p.attr({path:Raphael.transformPath(p.attr('path').toString(),'T'+x+','+y)});
		}
// if (p.id==12) console.log('move',shape,p);
	};
	var getChartPointCenter=function(p,shape='circle') {
		let x,y,r; // cx,cy, radius
		if (shape.match('circle|square|diamond')) {
			let midWidth;
			switch (shape) {
				case 'square':
					midWidth=p.attr('width')/2;
					x=p.attr('x')+midWidth;
					y=p.attr('y')+midWidth;
					r=midWidth;
					break;
				case 'diamond':
					midWidth=p.attr('width')/2;
					x=p.attr('x')+midWidth;
					y=p.attr('y')+midWidth;
					r=midWidth;
					break;
				default: // circle
					x=p.attr('cx');
					y=p.attr('cy');
					r=p.attr('r');
			}
		}
		else { // path
			x=p.attr('path')[0][1];
			y=p.attr('path')[0][2];
			r=Math.max(Math.abs(p.attr('path')[1][1]-x),Math.abs(p.attr('path')[1][2]-y)); // length of 1rst move: straight (up/down or right/left)
		}
// if (p.id==12) console.log('get',p);
		return [x,y,r];
	};


	/******************* Datasets visibility ***********************/
/*
	var setDatasetDisplay = function(mainC,dsIdx,visStatus) { // *public* used in HTML interface
		mainC.datasets[dsIdx].params.visible=visStatus;

		if (visStatus) {
			for (let i=0; i<mainC.datasets[dsIdx].data.length; i++) {
				let dp=mainC.datasets[dsIdx].data[i],
					C=(dp.subChart)? dp.subChart : mainC,
					minX=C.plotArea.attr('x'),
					maxX=minX+C.plotArea.attr('width')-1,
					minY=C.plotArea.attr('y'),
					maxY=minY+C.plotArea.attr('height')-1,
					p=dp.point || dp.svg;
				if (p.attr('cx') >= minX && p.attr('cx') <= maxX && p.attr('cy') >= minY && p.attr('cy') <= maxY) {
					p.show();
					//if (p.data('labelObj')) {p.data('labelObj').show()}
					if (p.data('showLabel')) {
						_emphasizePoint(p,'on',p.data('showLabel'));
					}
				}
			}
			if (mainC.datasets[dsIdx].dataConnection) {mainC.datasets[dsIdx].dataConnection.svg.show();} //connected points
		}
		else {
			for (let i=0; i<mainC.datasets[dsIdx].data.length; i++) {
				let dp=mainC.datasets[dsIdx].data[i];
				dp.point.hide();
				if (dp.point.data('labelObj')) {dp.point.data('labelObj').hide();}
				if (dp.pointList) { // multi-points
					for (let j=0; j<dp.pointList.length; j++) {
						dp.pointList[j].hide();
					}
				}
			}
			if (mainC.datasets[dsIdx].dataConnection) {mainC.datasets[dsIdx].dataConnection.svg.hide();} //connected points
		}
	};
*/	
	var setDatasetDisplay = function(mainC,chkbox,skipAutoExtend=false) { // *public* used in HTML interface
		let autoExtend=!skipAutoExtend; // default
		if (skipAutoExtend===false) {
			if (chkbox.readOnly) chkbox.checked=chkbox.readOnly=false;
			else if (!chkbox.checked) chkbox.readOnly=chkbox.indeterminate=true;
			let autoExtChk=document.getElementById('datasetAutoExtend_'+mainC.chartID);
			autoExtend=(autoExtChk && autoExtChk.checked)? true : false;
		}
		let visStatus=(chkbox.checked)? 2 : (chkbox.indeterminate)? 1 : 0,
			dsIdx=chkbox.value;
		mainC.datasets[dsIdx].params.visible=visStatus;

		/*Connection line*/
		let dsLine=mainC.datasets[dsIdx].dataConnection;
		if (dsLine) {
			if (visStatus==0) {dsLine.svg.hide();}
			else if (visStatus==1) {dsLine.svg.attr({stroke:idvis.dimColor,opacity:idvis.dimSetOpacity}).show();}
			else {dsLine.svg.attr({stroke:dsLine.color,opacity:dsLine.opacity}).show();}
//console.log(visStatus,dsLine.svg.attr('opacity'),dsLine);
		}
		
		/*Points*/
		const pointShape=mainC.datasets[dsIdx].params.pointShape;
		for (let i=0; i<mainC.datasets[dsIdx].data.length; i++) {
			let dp=mainC.datasets[dsIdx].data[i];
			if (!dp) continue; //categoryPlot with missing data
			let p=dp.point || dp.svg;
			if (visStatus==0) {p.hide();}
			else { // 1 or 2
				if (visStatus==1) {p.attr({fill:_findPointColor(dp),opacity:idvis.dimSetOpacity});}
				else {p.attr({fill:_findPointColor(dp),opacity:mainC.datasets[dsIdx].params.opacity});} //,stroke:C.datasets[dsIdx].params.color
				let C=(dp.subChart)? dp.subChart : mainC,
				    minX=C.plotArea.attr('x'),
					maxX=minX+C.plotArea.attr('width')-1,
					minY=C.plotArea.attr('y'),
					maxY=minY+C.plotArea.attr('height')-1,
					[pX,pY]=getChartPointCenter(p,pointShape);
				// if (p.attr('cx') >= minX && p.attr('cx') <= maxX && p.attr('cy') >= minY && p.attr('cy') <= maxY) { // }
				if (pX >= minX && pX <= maxX && pY >= minY && pY <= maxY) {
						p.show();
					//if (p.data('labelObj')) {p.data('labelObj').show()}
					if (visStatus==2 && p.data('showLabel')) {
						_emphasizePoint(p,'on',p.data('showLabel'),pX,pY);
					}
				}
			}
//if (i==0) console.log(visStatus,p.attr('opacity'),dp);
			if (visStatus <= 1) {
				//if (p.data('labelObj')) {p.data('labelObj').hide();} // unselected point
				if (p.data('showLabel')) {_selectPoint(p,'off');}
				if (dp.pointList) { // multi-points
					for (let j=0; j<dp.pointList.length; j++) {
						dp.pointList[j].hide();
					}
				}
			}
		}
		/* Auto-extend (selection status) to following checkboxes */
		if (autoExtend===true) {
			let dsChkBoxes=document.getElementsByName('datasetCheckboxes_'+mainC.chartID);
			let afterChkBox=false;
			for (let i=0; i<dsChkBoxes.length; i++) {
				if (dsChkBoxes[i].value==chkbox.value) {
					afterChkBox=true;
					continue;
				}
				else if (afterChkBox===false) {continue;}
				dsChkBoxes[i].checked=chkbox.checked;
				dsChkBoxes[i].readOnly=chkbox.readOnly;
				dsChkBoxes[i].indeterminate=chkbox.indeterminate;
				setDatasetDisplay(mainC,dsChkBoxes[i],true);
			}
		}
	};

	/******************* Point highlighting ***********************/
	var addHighlighting = function(mainC,hName,hColor,pointSet,matchPattern) { // *public* (use ### in matchPattern as point id)
		if (!mainC.allowHighlight) {
			alert('Highlighting is not allowed on chart #'+mainC.chartID);
			return false;
		}
		if (!mainC.highlightedPoints) {
			mainC.highlightedPoints={};
		}
		if (!mainC.highlightOrder) { // keeps track of order of layers of highlights
			mainC.highlightOrder=[];
		}
		if (mainC.highlightedPoints[hName]) {
			alert(hName+' is already used.');
			return false;
		}
		mainC.highlightedPoints[hName]={color:hColor,visible:true,pointSelected:false,dataPoints:[]};
		mainC.highlightOrder.push(hName);
		let dsIdxList={};
		for (let psIdx in pointSet) {
			if (psIdx==-1) { // use all datasets
				for (let i=0; i < mainC.datasets.length; i++) {dsIdxList[i]=-1;}
				break;
			}
			else {dsIdxList[psIdx]=psIdx;}
		}
		// Points update
		for (let dsIdx in dsIdxList) {
			if (!dsIdxList.hasOwnProperty(dsIdx)) {continue;}
			let psIdx=dsIdxList[dsIdx];
			P1:for (let i=0; i < pointSet[psIdx].length; i++) {
				let matchRegExp=(matchPattern)? new RegExp(matchPattern.replace('###',pointSet[psIdx][i])) : null;
				for (let j=0; j < mainC.datasets[dsIdx].data.length; j++) {
					if ((matchPattern && matchRegExp.exec(mainC.datasets[dsIdx].data[j].externalID)) || (!matchRegExp && mainC.datasets[dsIdx].data[j].externalID==pointSet[psIdx][i])) {
						mainC.highlightedPoints[hName].dataPoints.push(mainC.datasets[dsIdx].data[j]);
						mainC.datasets[dsIdx].data[j].highlightNames.push(hName);
						mainC.datasets[dsIdx].data[j].isVisble=true; // in case hideUnmatchedPoints is true & 1st active highlight
						let p=mainC.datasets[dsIdx].data[j].point;
						if (mainC.datasets[dsIdx].params.visible) {p.show();} // in case hideUnmatchedPoints is true & 1st active highlight
						p.toFront();
						if (p.data('showLabel')) {
							if (p.data('labelObj')) {
								_changeLabelColor(p.data('labelObj'),hColor);
								/*
								p.data('labelObj').remove();
								p.removeData('labelObj');
								_displayPointLabel(p,p.data('showLabel')); // label type
								*/
							}
						}
						else {p.attr('fill',hColor);}
						if (!matchPattern) continue P1; // more than 1 match per pointSet
					}
				}
			}
		}
		
		//Make sure all displayed labels above highlighted points
		_moveLabelsToFront(mainC);

		//Update legends on chart itself
		_updateHighlightLegends(mainC);

		//DIV update (after points because of point count display)
		_updateHighlightList(mainC);

		return true;
	};

	var deleteHighlighting = function(mainC,hName) { // *public*
		let hlPoints=mainC.highlightedPoints[hName].dataPoints;
		for (let i=0; i < hlPoints.length; i++) { // dataPoints
			//Removing this hl from list
			for (let j=0; j < hlPoints[i].highlightNames.length; j++) {
				if (hlPoints[i].highlightNames[j] === hName) {
					hlPoints[i].highlightNames.splice(j,1);
				}
			}
			let newHlNames=hlPoints[i].highlightNames;
			let p=hlPoints[i].point;
			if (mainC.hideUnmatchedPoints && newHlNames.length===0) {
				hlPoints[i].isHidden=true;
				p.hide();
			}
			if (p.data('showLabel')) { // point is selected: color=#f00. Do not change it!
				if (p.data('labelObj')) {
					_changeLabelColor(p.data('labelObj'),_findPointColor(hlPoints[i],true));
					/*
					p.data('labelObj').remove();
					p.removeData('labelObj');
					if (!hlPoints[i].isHidden) {_displayPointLabel(p,p.data('showLabel'));} // change only label color
					*/
				}
			}
			else {
				p.attr('fill',_findPointColor(hlPoints[i]));
			}
		}
		delete mainC.highlightedPoints[hName];
		let mustRefresh=(mainC.highlightOrder[0]===hName)? false : true;
		let idx=mainC.highlightOrder.findIndex(hl => hl===hName);
		mainC.highlightOrder.splice(idx,1); // remove hName from list
		
		if (mainC.hideUnmatchedPoints) {
			let noHighlighting=true;
			for (let hlName in mainC.highlightedPoints) {
				if (mainC.highlightedPoints.hasOwnProperty(hlName)) {noHighlighting=false; break;}
			}
			if (noHighlighting) { // reset hideUnmatchedPoints to false
				for (let dsIdx=0; dsIdx < mainC.datasets.length; dsIdx++) {
					for (let i=0; i < mainC.datasets[dsIdx].data.length; i++) {
						mainC.datasets[dsIdx].data[i].isHidden=false;
						if (mainC.datasets[dsIdx].params.visible) mainC.datasets[dsIdx].data[i].point.show();
					}
				}
			}
			mainC.hideUnmatchedPoints=false;
		}
		
		//Redraw HL to move de-highlighted to back
		if (mustRefresh) _redrawHighlightings(mainC);

		//Update legends on chart itself
		_updateHighlightLegends(mainC);

		//DIV update
		_updateHighlightList(mainC);

		//Callback function
		if (mainC.updateHighlight) {
			mainC.updateHighlight.callback(hName,'delete');
		}
	};

	var editHighlighting = function(mainC,action,hlIndex) { // *public*
		if (action==='edit') {
			document.getElementById(mainC.divID+'_hlightDispSPAN'+hlIndex).style.display='none';
			document.getElementById(mainC.divID+'_hlightEditSPAN'+hlIndex).style.display='';
		}
		else {
			document.getElementById(mainC.divID+'_hlightDispSPAN'+hlIndex).style.display='';
			document.getElementById(mainC.divID+'_hlightEditSPAN'+hlIndex).style.display='none';
		}
	};

	var applyHighlightEdition = function(mainC,hlIndex) { // *public*
		let hlInpName=document.getElementById(mainC.divID+'_hlightNewName'+hlIndex),
			newName=hlInpName.value,
			oldName=hlInpName.dataset.oldName;
		if (newName !== oldName) { // check for duplicate
			for (let hlName in mainC.highlightedPoints) {
				if (newName===oldName) {
					alert(newName+' is already used.');
					return;
				}
			}
		}
		let oldColor=mainC.highlightedPoints[oldName].color.toUpperCase(),
			newColor=document.getElementById(mainC.divID+'_hlightNewColor'+hlIndex).value.toUpperCase(),
			hasChanged=true; // default
		if (newName===oldName) {
			if (newColor===oldColor) { // nothing to update
				//alert('No change was made to '+oldName);
				//return; <--- No return to allow update even if no change due to Firefox loosing control of input color if mulitple hl edited at once
				hasChanged=false;
			}	
			else {mainC.highlightedPoints[oldName].color=newColor;}
		}
		else {
			mainC.highlightedPoints[newName]={color:newColor,visible:mainC.highlightedPoints[oldName].visible,dataPoints:[]};
			let hlPoints=mainC.highlightedPoints[oldName].dataPoints;
			for (let i=0; i < hlPoints.length; i++) { // dataPoints
				for (let j=0; j < hlPoints[i].highlightNames.length; j++) {
					if (hlPoints[i].highlightNames[j]===oldName) { //Replace oldName with newName
						hlPoints[i].highlightNames[j]=newName;
						break;
					}
				}
				mainC.highlightedPoints[newName].dataPoints.push(hlPoints[i]); // add data point to new
			}
			mainC.highlightedPoints[oldName].dataPoints=[];
			delete mainC.highlightedPoints[oldName];
			let idx=mainC.highlightOrder.findIndex(hl => hl===oldName);
			mainC.highlightOrder[idx]=newName;
		}

		// Update points color
		if (newColor !== oldColor && document.getElementById(mainC.divID+'_hlightCHK'+hlIndex).checked) { // hl is active
			setHighlightVisibility(mainC,newName,true);
		}

        //Update legends on chart itself
        _updateHighlightLegends(mainC);

        //DIV update (after points because of point count display)
        _updateHighlightList(mainC);

        if (mainC.updateHighlight && mainC.updateHighlight.callback) {mainC.updateHighlight.callback(oldName,'edit',newName,newColor);}
	};
	
	var selectHighlightedPoints = function(mainC,hName) { // *public*
		let hl=mainC.highlightedPoints[hName],
			action=(hl.pointSelected===true)? 'off' : 'on';
		for (let i=0; i<hl.dataPoints.length; i++) {
			if (!hl.dataPoints[i].isHidden) {
				_selectPoint(hl.dataPoints[i].point,action,'min');
			}
		}
		hl.pointSelected=!hl.pointSelected;
		if (hl.pointSelected) _moveLabelsToFront(mainC);
    };

    var setHighlightVisibility = function(mainC,hName,visStatus) { // *public*
        mainC.highlightedPoints[hName].visible=visStatus;
        let hlPoints=mainC.highlightedPoints[hName].dataPoints;
        for (let i=0; i < hlPoints.length; i++) { // dataPoints
            let dp=hlPoints[i];
            //if (dp.highlightNames[dp.highlightNames.length-1] != hName) {continue;} // no effects on point color
            if (visStatus===false && !dp.isHidden && mainC.hideUnmatchedPoints) {
                let hidePoint=true;
                for (let j=0; j < dp.highlightNames.length; j++) {
                    if (dp.highlightNames[j] != hName && mainC.highlightedPoints[dp.highlightNames[j]].visible) {
                        hidePoint=false;
                        break;
                    }
                }
                dp.isHidden=hidePoint;
            }
            let p=dp.point;
            if (p.data('showLabel')) { // point is selected: color=#f00. Do not change it!
                if (dp.isHidden) {
                    _selectPoint(p,'off');
                    //p.hide();
                }
                else if (p.data('labelObj')) {
					_changeLabelColor(p.data('labelObj'),_findPointColor(dp,true));
					/*
                    p.data('labelObj').remove();
                    p.removeData('labelObj');
                    _displayPointLabel(p,p.data('showLabel')); // change only label color
                    //if (!dp.isHidden) _displayPointLabel(p,p.data('showLabel')); // change only label color
					*/
                }
            }
            else {
                // highlight is the newest for this point (&& other highlights exist) => update point color
                //let newColor;
                if (visStatus===true) { // show highlight
                    for (let h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
                        if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
                            //newColor=mainC.highlightedPoints[dp.highlightNames[h]].color;
                            dp.isHidden=false;
                            break;
                        }
                    }
                    //p.show();
                }
				/*
                else { // hide highlight
                    newColor=dp.dataset.params.color; // default
                    for (let h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
                        if (dp.highlightNames[h]==hName) {continue;}
                        if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
                            newColor=mainC.highlightedPoints[dp.highlightNames[h]].color;
                            break;
                        }
                    }
                    //if (dp.isHidden) p.hide();
                }
                */
                p.attr('fill',_findPointColor(dp));
            }
			if (dp.dataset.params.visible) {
				if (dp.isHidden) {p.hide();} else {p.show();}
			}
			if (visStatus===true) p.toFront();
        }

		//Make sure all displayed labels above highlighted points
		if (visStatus === true) {  // move hName above all others
			for (let i=0; i<mainC.highlightOrder.length; i++) {
				if (mainC.highlightOrder[i]===hName) {
					mainC.highlightOrder.splice(i,1); // remove hName from list
					break;
				}
			}
			mainC.highlightOrder.push(hName);
			_moveLabelsToFront(mainC);
		}
		else { // move other highlightings above these "unhighlighted" points
			let lastHiddenIdx=-1;
			for (let i=0; i<mainC.highlightOrder.length; i++) {
				if (mainC.highlightOrder[i]===hName) {
					mainC.highlightOrder.splice(i,1); // remove hName from list (reinserted below)
					break;
				}
				else if (!mainC.highlightedPoints[ mainC.highlightOrder[i] ].visible) {
					lastHiddenIdx=i; // records latest non-visible
				}
			}
			mainC.highlightOrder.splice(lastHiddenIdx+1,0,hName); // insert after last non-visible
			_redrawHighlightings(mainC);
		}

		//Update legends on chart itself
		_updateHighlightLegends(mainC);
	};

	var setUnmatchedPointsVisibility = function(mainC,hideStatus) { // *public*
		mainC.hideUnmatchedPoints=hideStatus;
		for (let dsIdx=0; dsIdx < mainC.datasets.length; dsIdx++) {
			POINT:for (let i=0; i < mainC.datasets[dsIdx].data.length; i++) {
				let dp=mainC.datasets[dsIdx].data[i];
				if (hideStatus===true) { // hide unmatched points
					for (let h=0; h < dp.highlightNames.length; h++) {
						if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
							continue POINT;
						}
					}
					dp.isHidden=true;
					if (dp.point.data('showLabel')) _selectPoint(dp.point,'off');
					dp.point.hide();
				}
				else {
					dp.isHidden=false;
					if (mainC.datasets[dsIdx].params.visible) dp.point.show();
				}
			}
		}
	};

	function _redrawHighlightings(mainC) {
		mainC.highlightOrder.forEach(function(hName) {
			if (!mainC.highlightedPoints[hName].visible) {return;} // (<=> "continue" in for loop) hidden highlighting
			let hlPoints=mainC.highlightedPoints[hName].dataPoints;
			for (let i=0; i < hlPoints.length; i++) {
				hlPoints[i].point.toFront();
			}
		});
		_moveLabelsToFront(mainC);
	}
	
    function _updateHighlightLegends(mainC) { // private
		let paper=mainC.canvas;
		if (!mainC.legendX) {
			mainC.highlightLegends=paper.set();  // initialize SVG set
			mainC.legendX=paper.width; // record initial width
		}
		else {mainC.highlightLegends.remove();} // clear previous legend if any
		let posY=25,
			numAnnot=0;
		mainC.highlightOrder.forEach(function(hName) { // mainC.highlightedPoints
			if (!mainC.highlightedPoints[hName].visible) {return;} // (<=> "continue" in for loop) hidden highlighting
			numAnnot++;
			if (numAnnot==1) {
				mainC.highlightLegends.push(paper.text(mainC.legendX,posY,'Legends:').attr({'font-size':12,'font-weight':'bold','text-anchor':'start'}));
				posY+=17;
			}
			mainC.highlightLegends.push(paper.rect(mainC.legendX+5,posY-5,10,10).attr({fill:mainC.highlightedPoints[hName].color,stroke:mainC.highlightedPoints[hName].color}));
			mainC.highlightLegends.push(paper.text(mainC.legendX+20,posY,hName).attr({'font-size':12,'text-anchor':'start'}));
			posY+=15;
		});
		/*Readjust chart size if necessary */
		let newWidth=Math.max(mainC.legendX,mainC.highlightLegends.getBBox().x2+15),
			oldWidth=paper.width;
		paper.setSize(newWidth,paper.height);
		paper.bottom.attr({width:newWidth}); // adjust background panel too
		if (oldWidth > newWidth) { // make sure that all labels are still visible
			var redrawn=false;
			for (let dsIdx=0; dsIdx<mainC.datasets.length; dsIdx++) {
				if (mainC.datasets[dsIdx].params.visible < 2) continue; // must be fully visible
				for (let i=0; i < mainC.datasets[dsIdx].data.length; i++) {
					let dp=mainC.datasets[dsIdx].data[i];
					let p=dp.point;
					if (p.data('labelObj') && !dp.isHidden) {
						let tBox=p.data('labelObj')[0][0];
						if (tBox.attr('x')+tBox.attr('width') > newWidth) { // label is not fully visible +> redraw
							p.data('labelObj').remove();
							p.removeData('labelObj');
							_displayPointLabel(p,p.data('showLabel')); // label type
							redrawn=true;
						}
					}
				}
			}
			if (redrawn===true) {
				_moveDragAreasToFront(mainC);
			}
		}
    }

    function _updateHighlightList(mainC) { // private DIV content update
		// let hlNameList=mainC.highlightOrder.sort(); // list in form is ordered alphabetically
        let hlNameList=[...mainC.highlightOrder]; // Copy!
        hlNameList.sort(); // list in form is ordered alphabetically
        let hCode='';
        if (hlNameList.length) {
            hCode='<INPUT type="checkbox" value="1" onclick="idvis.lib.setUnmatchedPointsVisibility(idvis.registeredCharts['+mainC.chartID+'],this.checked)"';
            if (mainC.hideUnmatchedPoints) {hCode+=' checked';}
            hCode+='>Hide unmatched points<BR>';
			const buttonStyleStrg='border:none;background-color:Transparent;font-size:14px;font-weight:bold;padding:0px;margin:0px 2px';
            for (let i=0; i<hlNameList.length; i++) {
                let hlName=hlNameList[i];
				if (mainC.updateHighlight && mainC.updateHighlight.editable) { // SPAN covers hl & buttons
					hCode+='<SPAN id="'+mainC.divID+'_hlightDispSPAN'+i+'">';
				}
                hCode+='<INPUT type="checkbox" id="'+mainC.divID+'_hlightCHK'+i+'" value="'+hlName+'" onclick="idvis.lib.setHighlightVisibility(idvis.registeredCharts['+mainC.chartID+'],this.value,this.checked)"';
                if (mainC.highlightedPoints[hlName].visible) hCode+=' checked';
                hCode+='><A href="javascript:idvis.lib.selectHighlightedPoints(idvis.registeredCharts['+mainC.chartID+'],\''+hlName+'\')" style="font-size:12px;color:'+mainC.highlightedPoints[hlName].color+'">'+hlName+' ('+mainC.highlightedPoints[hlName].dataPoints.length+')</A>&nbsp;';
                if (!mainC.updateHighlight || !mainC.updateHighlight.editable) { // SPAN covers only edit & delete buttons
					hCode+='<SPAN id="'+mainC.divID+'_hlightDispSPAN'+i+'">';
				}
				hCode+='<INPUT type="button" value="&rarr;" style="'+buttonStyleStrg+'" onclick="idvis.lib.editHighlighting(idvis.registeredCharts['+mainC.chartID+'],\'edit\','+i+')">'; // \''+hlName+'\',
				//hCode+='<BUTTON type="button" style="padding:0px;width:20px" onclick="idvis.lib.editHighlighting(idvis.registeredCharts['+mainC.chartID+'],\'edit\',\''+hlName+'\')"><img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAYAAAAfSC3RAAAAXUlEQVQokc3QQQ3AIAyF4V8CEpAyKUiYJJwMB8MRO/CaHddyWHhJ0/TwpU3hhxxRkIAbGOphZPW52dAJdKHuRUWz4RRBZS8E0FZQ5v2YG8H82qVTmxcBVOHsBct5APW0HZKx7Nn3AAAAAElFTkSuQmCC"></BUTTON>';
				// &larr; = <- , &otimes; = x in circle , &check; , &cross; &olcross; &plusmn;
                hCode+='<INPUT type="button" value="&cross;" style="'+buttonStyleStrg+'" onclick="idvis.lib.deleteHighlighting(idvis.registeredCharts['+mainC.chartID+'],\''+hlName+'\')"><BR>';
				hCode+='</SPAN>';
				hCode+='<SPAN id="'+mainC.divID+'_hlightEditSPAN'+i+'" style="display:none">';
				if (mainC.updateHighlight && mainC.updateHighlight.editable) {
					hCode+='<INPUT type="text" id="'+mainC.divID+'_hlightNewName'+i+'" value="'+hlName+'" data-old-name="'+hlName+'" style="width:120px">'; // data-old-name -> .dataset.oldName
				}
				else { // Do not allow hlName to be edited
					hCode+='<INPUT type="hidden" id="'+mainC.divID+'_hlightNewName'+i+'" value="'+hlName+'" data-old-name="'+hlName+'">'; // data-old-name -> .dataset.oldName
				}
				hCode+='<INPUT type="color" name="'+mainC.divID+'_hlightNewColor'+i+'" id="'+mainC.divID+'_hlightNewColor'+i+'" value="'+mainC.highlightedPoints[hlName].color+'" style="height:20px;width:50px;vertical-align:middle">';
				hCode+='<INPUT type="button" value="&check;" style="'+buttonStyleStrg+'" onclick="idvis.lib.applyHighlightEdition(idvis.registeredCharts['+mainC.chartID+'],'+i+')"><INPUT type="button" value="&larr;" style="'+buttonStyleStrg+'" onclick="idvis.lib.editHighlighting(idvis.registeredCharts['+mainC.chartID+'],\'cancel\','+i+')">';
				hCode+='</SPAN>';
			}
        }
        else {hCode='None';}

        document.getElementById(mainC.divID+'_hlight').innerHTML=hCode;
    }


    /******************* Search ***********************/
	var searchJobs={};
    var searchDataPoints = function(mainC,searchText,searchResDivID,extSearchChkID) { // *public*
		var searchResDiv=document.getElementById(searchResDivID);
        if (!searchText || searchText.length < 2) {
            searchResDiv.innerHTML='<FONT style="color:#DD0000">Search string is too short!</FONT>';
            return;
        }
		if (extSearchChkID !== null && document.getElementById(extSearchChkID).checked) {
			mainC.extSearchJobID=Date.now();
			searchJobs[mainC.extSearchJobID]=mainC;
			searchResDiv.innerHTML='<FONT style="color:#0000FF">Processing...</FONT>';
			mainC.searchable.externalSearch(searchText,mainC.extSearchJobID,processExternalSearchResult); // (local search function,array of function arguments,index of search text in array). To be recalled by external search function
			return;
		}
        let matchList=[],
			matchExp=new RegExp(searchText,"i"),
			neverMatched=true, okProceed100=false, okProceed1000=false,
			newZoom={};
        for (let dsIdx=0; dsIdx<mainC.datasets.length; dsIdx++) {
			if (!mainC.datasets[dsIdx].params.visible) continue;
			const pointShape=mainC.datasets[dsIdx].params.pointShape;
            for (let i=0; i < mainC.datasets[dsIdx].data.length; i++) {
                var dp=mainC.datasets[dsIdx].data[i];
                if (dp.isHidden) continue;
                if (dp.label.match(matchExp)) {
                    matchList.push(dp);
                    if (matchList.length > 100 && !okProceed100) {
                        if (!confirm('More than 100 matches found! Proceed?')) {
                            searchResDiv.innerHTML='<FONT style="color:#DD0000">More than 100 matches found!</FONT>';
                            return;
                        }
                        okProceed100=true;
                    }
                    if (matchList.length > 1000 && !okProceed1000) {
                        if (!confirm('More than 1000 matches found! Proceed?')) {
                            searchResDiv.innerHTML='<FONT style="color:#DD0000">More than 1000 matches found!</FONT>';
                            return;
                        }
                        okProceed1000=true;
					}
					let [pX]=getChartPointCenter(dp.point,pointShape);
					// if (dp.point.attr('cx') < 0) { //} point if off chart ZoomOut required
					if (pX < 0) { // point if off chart ZoomOut required
                        let chartName=(dp.subChart)? dp.subChart.name : 'main';
                        newZoom[chartName]=true;
                    }
                    neverMatched=false;
                }
            }
        }
        if (neverMatched) {
            //alert('No match found!');
            searchResDiv.innerHTML='<FONT style="color:#DD0000">No match found!</FONT>';
            return;
        }
        for (let chName in newZoom) {
            if (chName=='main') {zoomOut(mainC,true);}
            else {zoomOut(mainC.subChart[chName],true);}
        }
        for (let i=0; i < matchList.length; i++) {
            //_emphasizePoint(matchList[i].point,'on','min');
            _selectPoint(matchList[i].point,'on','min');
        }
		_moveLabelsToFront(mainC); // make sure all displayed labels stay above matched points
        searchResDiv.innerHTML=matchList.length+' match(es) found!';
    };

	var processExternalSearchResult = function(extSearchJobID,searchText) {
		var mainC=searchJobs[extSearchJobID];
		if (!mainC || !mainC.extSearchJobID || mainC.extSearchJobID != extSearchJobID) {return;} // ignore this job (another search was launched after it return)
		delete searchJobs[mainC.extSearchJobID]; // ext job is done
		mainC.extSearchJobID=undefined; // ext job is done
		var searchResDivID=mainC.divID+'_srchRes';
		searchDataPoints(mainC,searchText,searchResDivID,null);
	};
	
    /******************* Selection management ***********************/
    var manageSelection = function(mainC,action) { // *public*
        let chartList=[];
        if (mainC.subChart) {
            for (let c=0; c<mainC.activeCharts.length; c++) {
                chartList.push(mainC.subChart[mainC.activeCharts[c]]);
            }
        }
        else {chartList.push(mainC);}

        // Looping through displayed charts
        for (let c=0; c<chartList.length; c++) {
            var C=chartList[c];
            if (!C.dragArea || C.dragArea.data('status')=='off') continue; //return;
            let minPx=C.dragArea.attr('x'),
				maxPx=minPx+C.dragArea.attr('width')-1,
				minPy=C.dragArea.attr('y'),
				maxPy=minPy+C.dragArea.attr('height')-1;
            for (let dsIdx=0; dsIdx<mainC.datasets.length; dsIdx++) {
				if (!mainC.datasets[dsIdx].params.visible) continue;
				const pointShape=mainC.datasets[dsIdx].params.pointShape;
                for (let i=0; i < mainC.datasets[dsIdx].data.length; i++) {
                    let dp=mainC.datasets[dsIdx].data[i];
                    if (dp.isHidden || (dp.subChart && dp.subChart !== C)) continue;
					let p=dp.point;
					let [pX,pY]=getChartPointCenter(p,pointShape);
					// if (p.attr('cx')>=minPx && p.attr('cx')<=maxPx && p.attr('cy')>=minPy && p.attr('cy')<=maxPy) { //}
					if (pX>=minPx && pX<=maxPx && pY>=minPy && pY<=maxPy) {
                        _selectPoint(p,action,'min',pX,pY);
                        //_emphasizePoint(p,action,'min');
                    }
                }
            }
        }
    };

    var listSelectedPoints = function(C,actionIndex) { // *public*
        let selectedPoints={},
			count=0;
        for (let dsIdx=0; dsIdx<C.datasets.length; dsIdx++) {
			if (!C.datasets[dsIdx].params.visible) continue;
			const pointShape=C.datasets[dsIdx].params.pointShape;
            selectedPoints[dsIdx]=[];
            for (let i=0; i < C.datasets[dsIdx].data.length; i++) {
				let dp=C.datasets[dsIdx].data[i];
				let [pX]=getChartPointCenter(dp.point,pointShape);
				// if (dp.point.attr('cx') > 0 && dp.point.data('labelObj')) { //} visible & selected
				if (pX > 0 && dp.point.data('labelObj')) { // visible & selected
                    selectedPoints[dsIdx].push(dp.externalID);
                    count++;
                }
            }
            if (selectedPoints[dsIdx].length===0) {delete selectedPoints[dsIdx];}
        }
        if (count===0) {alert('No selected points found in displayed range.');}
        else {
            let thresholds={};
            for (let i=0; i<C.features.length;i++) { // also export list of thresholds name:value
                if (C.features.type==='line' && C.features[i].name) {
					thresholds[C.features[i].name]=C.features[i].initValue;
				}
            }
            if (actionIndex < 0) {
                C.pointOnList(selectedPoints,thresholds);
            }
            else {
                C.pointOnList[actionIndex][1](selectedPoints,thresholds);
            }
        }
    };

    /******************* Thresholds management ***********************/
    var selectThreshold = function(C,thIdx,thresBoxID) { // *public*
        document.getElementById(thresBoxID).value=(thIdx>=0)? C.editableThresholds[thIdx].initValue : '';
    };

    var updateThreshold = function(mainC,thresSelID,thresBoxID) { // *public*
        let thIdx=document.getElementById(thresSelID).value,
			newValue=document.getElementById(thresBoxID).value;
        if (!idvis.isNumber(newValue)) {alert(newValue+' is not a valid number!'); return;}
        if (thIdx<0 || !newValue) return;
		let thLine=mainC.editableThresholds[thIdx];
        thLine.setValues(newValue);
        //thLine.path.hide(); // No animation
        let chartX=thLine.chart.plotArea.attr('x'), chartX0=chartX-0.5,
			chartY=thLine.chart.plotArea.attr('y'), chartY0=chartY-1.5,
			chartW=thLine.chart.plotArea.attr('width'),
			chartH=thLine.chart.plotArea.attr('height'),
			settings=thLine.chart.chartSettings.current,
			prevPos=thLine.pathPos,
			axisMin,axisMax;
		if (mainC.chartType==='categoryPlot') { // TODO: Find a more universal solution
			chartX0+=0.5;
			chartY0+=1.5;
		}
		if (thLine.axis.match('X')) {
			axisMin=chartX; axisMax=chartX0+chartW;
			thLine.pathPos=chartX0+Math.round((thLine.value-settings[thLine.axis].minRange)/settings[thLine.axis].pix2valRatio);
		}
		else {
			axisMin=chartY; axisMax=chartY0+chartH;
			if (mainC.horizontalCP) { // Horizontal CateogryPlot Y0 is at top of chart
				thLine.pathPos=axisMin+Math.round((thLine.value-settings[thLine.axis].minRange)/settings[thLine.axis].pix2valRatio);
			}
			else { // usual cases
				thLine.pathPos=axisMax-Math.round((thLine.value-settings[thLine.axis].minRange)/settings[thLine.axis].pix2valRatio);
			}
		}
		/* No animation */
		//if (thLine.pathPos >= axisMin && thLine.pathPos <= axisMax) {
		//	if (thLine.axis.match('X')) {thLine.path.translate(thLine.pathPos-prevPos,0);}
		//	else {thLine.path.translate(0,thLine.pathPos-prevPos);}
		//	thLine.path.show();
		//}
		/* Animation version */
		let m = Raphael.matrix(1, 0, 0, 1, 0, 0); /* the identity matrix */
		if (thLine.axis.match('X')) {m.translate(thLine.pathPos-prevPos,0);}
		else {m.translate(0,thLine.pathPos-prevPos);}
		if (thLine.pathPos >= axisMin && thLine.pathPos <= axisMax) {thLine.path.show();}
		else {thLine.path.hide();}
		let newPath = Raphael.mapPath(thLine.path.attr("path"), m);
		thLine.path.animate({"path": newPath},200,"easeInOut");
		
		/* Chart-Specific callback */
		if (mainC.onThresholdUpdate) {mainC.onThresholdUpdate();}
    };
	
	var activateDimArea = function(mainC,status) {
		mainC.dimArea.active=status;
		/* Chart-Specific callback */
		if (mainC.onThresholdUpdate) {mainC.onThresholdUpdate();}
	};

	/******************* Color selection ***********************/
	var toggleColorSelection = function(mainC,elementStrg) {
		var colorDiv=document.getElementById('colorDIV_'+mainC.divID);
		if (elementStrg==='none') { // nothing selected
			colorDiv.style.display='none';
			return;
		}
		colorDiv.style.display='';
		const [elemType,idx]=elementStrg.split(':');
		var curColor=(elemType==='dataset')? mainC.datasets[idx].params.color : mainC.features[idx].color;
		//if (elemType==='dataset') {curColor=mainC.datasets[idx].params.color;}
		//else if (elemType==='feature') {curColor=mainC.features[idx].color;}
		if (curColor.length <= 4) { // convert #ABC to #AABBCC
			curColor=curColor.replace('#','');
			curColor = curColor.split('').map(function (hex) {return hex + hex;}).join('');
			curColor='#'+curColor;
		}
		document.getElementById('colorPicker_'+mainC.chartID).value=curColor;
	};

	var changeElementColor = function(mainC,elementStrg,newColor) {
//console.log(mainC.divID,elementStrg,newColor);
		const [elemType,idx]=elementStrg.split(':');
		var curColor=(elemType==='dataset')? mainC.datasets[idx].params.color : mainC.features[idx].color;
		if (newColor==curColor) return;
		/* Dataset */
		if (elemType==='dataset') {
			mainC.datasets[idx].params.color=newColor;
			//Icon
			if (mainC.datasets[idx].params.icon) {mainC.datasets[idx].params.icon.attr({fill:newColor});}
			//Name displayed
			var dsName=document.getElementById('dataset_'+mainC.chartID+'_'+idx);
			if (dsName) dsName.style.color=newColor;

			//Points
			let chkbox=document.getElementById('datasetChkbox_'+mainC.chartID+'_'+idx),
				visStatus=(!chkbox || chkbox.checked)? 2 : (chkbox.indeterminate)? 1 : 0;
			if (visStatus < 2) return; // dim state: nothing to do

			for (let i=0; i<mainC.datasets[idx].data.length; i++) {
				let dp=mainC.datasets[idx].data[i];
				if (!dp) continue; //categoryPlot with missing data
				let p=dp.point || dp.svg,
					pColor=_findPointColor(dp);
//console.log('dp',i,pColor,p.attr('fill'));
				if (pColor==p.attr('fill')) continue; // point does not change color
				p.attr({fill:pColor});
				//Point label
				if (p.data('showLabel')) {
					_emphasizePoint(p,'on',p.data('showLabel'));
					if (p.data('labelObj')) {
						_changeLabelColor(p.data('labelObj'),pColor);
					}
				}
			}
		}
		else if (elemType==='feature') {
			mainC.features[idx].color=newColor;
			let attribute=(mainC.features[idx].type==='text')? 'fill' : 'stroke';
			mainC.features[idx].path.attr(attribute,newColor);
			//...
		}


	};
	
    /******************* Chart zooming ***********************/
    var zoomIn = function(C) { // *public* (not actually used as public, but could be one day)
        //let mainC=(C.mainChart)? C.mainChart : C;
        let chartX0=C.plotArea.attr('x')-1,
			chartY0=C.plotArea.attr('y')-1,
			chartW=C.plotArea.attr('width'),
			chartH=C.plotArea.attr('height'),
			minPx,maxPx,minPy,maxPy;
        if (C.noZoomX) {
            minPx=chartX0;
            maxPx=chartX0+chartW;
        }
        else {
            minPx=C.dragArea.attr('x');
            maxPx=minPx+C.dragArea.attr('width')-1;
        }
        if (C.noZoomY) {
            minPy=chartY0;
            maxPy=chartY0+chartH;
        }
        else {
            minPy=C.dragArea.attr('y');
            maxPy=minPy+C.dragArea.attr('height')-1;
        }
        let ref=C.chartSettings.reference;
        let cur=C.chartSettings.current;

        for (let axis in C.usedAxes) {
			if (!C.usedAxes.hasOwnProperty(axis)) {continue;}
            let minV,maxV,axisSize;
            if (axis.match('X')) { // X,X2
				if (C.noZoomX) {continue;}
                minV=cur[axis].minRange+((minPx-chartX0)*cur[axis].pix2valRatio);
                maxV=cur[axis].minRange+((maxPx-chartX0)*cur[axis].pix2valRatio);
                axisSize=chartW;
            }
            else { // Y,Y2
				if (C.noZoomY) {continue;}
                minV=cur[axis].minRange+((chartY0+chartH-maxPy)*cur[axis].pix2valRatio); // maxPy for minVy !!!
                maxV=cur[axis].minRange+((chartY0+chartH-minPy)*cur[axis].pix2valRatio);
                axisSize=chartH;
            }
            let newRange=getChartScaleRange(false,minV,maxV,axisSize);

            cur[axis].minRange=minV;
            cur[axis].optMinRange=newRange[1];
            cur[axis].maxRange=maxV;
            cur[axis].tickSize=newRange[3];
            cur[axis].pix2valRatio=(cur[axis].maxRange-cur[axis].minRange)/axisSize;

            let axisZoom=(axis.match('X'))? 'zoomX' : 'zoomY';
            cur[axisZoom]=ref[axis].pix2valRatio/cur[axis].pix2valRatio;
        }

        /*** Clear & plot chart ***/
        clearChart(C);
        plotChart(C);
    };

    var zoomOut = function(C,zoomTo1) { // *public*
        let cur=C.chartSettings.current;
        if (cur.zoomX==1 && cur.zoomY==1) {return;}
        let ref=C.chartSettings.reference;

        for (let axis in C.usedAxes) {
			if (!C.usedAxes.hasOwnProperty(axis)) {continue;}
            let minV,maxV,axisSize,axisZoom;
            if (idvis.modKeyPressed || zoomTo1) { // set zooms to 1
                minV=ref[axis].minRange;
                maxV=ref[axis].maxRange;
            }
            else {
                minV=Math.max(cur[axis].minRange-(cur[axis].maxRange-cur[axis].minRange)/2,ref[axis].minRange);
                maxV=Math.min(cur[axis].maxRange+(cur[axis].maxRange-cur[axis].minRange)/2,ref[axis].maxRange);
            }
            if (axis.match('X')) {
                if (C.noZoomX) {continue;}
                axisSize=C.plotArea.attr('width');
                axisZoom='zoomX';
            }
            else {
                if (C.noZoomY) {continue;}
                axisSize=C.plotArea.attr('height');
                axisZoom='zoomY';
            }
            if (minV==ref[axis].minRange && maxV==ref[axis].maxRange) { // zoomX=1
                cur[axis].minRange=ref[axis].minRange;
                cur[axis].optMinRange=ref[axis].optMinRange;
                cur[axis].maxRange=ref[axis].maxRange;
                cur[axis].tickSize=ref[axis].tickSize;
                cur[axis].pix2valRatio=ref[axis].pix2valRatio;
                cur[axisZoom]=1;
            }
            else {
                var newRange=getChartScaleRange(false,minV,maxV,axisSize);
                cur[axis].minRange=minV;
                cur[axis].optMinRange=newRange[1];
                cur[axis].maxRange=maxV;
                cur[axis].tickSize=newRange[3];
                cur[axis].pix2valRatio=(cur[axis].maxRange-cur[axis].minRange)/axisSize;
                cur[axisZoom]=ref[axis].pix2valRatio/cur[axis].pix2valRatio;
            }
        }

        /*** Clear & plot chart ***/
        clearChart(C);
        plotChart(C);
    };


    /******************* Chart panning ***********************/
    var panChart = function(C) { // *public*
        let chartW=C.plotArea.attr('width'),
            chartH=C.plotArea.attr('height'),
            ref=C.chartSettings.reference,
            cur=C.chartSettings.current;

        for (let axis in C.usedAxes) {
			if (!C.usedAxes.hasOwnProperty(axis)) {continue;}
            let dV,axisSize;
            if (axis.match('X')) { // opposite movement
                if (C.noZoomX) {continue;}
                dV=-C.panProcess.dx*cur[axis].pix2valRatio;
                axisSize=chartW;
            }
            else { // counter-opposite movement
                if (C.noZoomY) {continue;}
                dV=C.panProcess.dy*cur[axis].pix2valRatio;
                axisSize=chartH;
            }
            dV=(dV > 0)? Math.min(dV,ref[axis].maxRange-cur[axis].maxRange) : Math.max(dV,ref[axis].minRange-cur[axis].minRange);
            cur[axis].minRange+=dV;
            cur[axis].maxRange+=dV;
            let newRange=getChartScaleRange(false,cur[axis].minRange,cur[axis].maxRange,axisSize);
            cur[axis].optMinRange=newRange[1];
            cur[axis].tickSize=newRange[3];
        }

        /*** Clear & plot chart ***/
        clearChart(C);
        plotChart(C);
    };

    var clearChart = function(C) { // *public*
        _clearDragArea(C.dragArea);

        if (C.features) {
            for (let i=0; i<C.features.length; i++) {
                if (C.features[i].path) {
                    C.features[i].path.hide();
                    /**/
                    C.features[i].path.remove();
                    C.features[i].path=null;
                    C.features[i].pathPos=null;
                    if (C.features[i].popup) {
                        C.features[i].popup.remove();
                        C.features[i].popup=null;
                    }
                    /**/
                }
            }
        }
        for (let i=0; i<C.chartMarks.length; i++) {
            C.chartMarks[i].remove();
        }
        C.chartMarks.length=0;

    //var maxDx=C.plotArea.attr('width')+50;
    //var maxDy=C.plotArea.attr('height')+50;

        let mainC=(C.mainChart)? C.mainChart : C;
        if (mainC.datasets) {
            for (let dsIdx=0; dsIdx<mainC.datasets.length; dsIdx++) {
				const pointShape=mainC.datasets[dsIdx].params.pointShape;
                for (let i=0; i<mainC.datasets[dsIdx].data.length; i++) {
                    let dp=mainC.datasets[dsIdx].data[i];
                    if (dp.subChart && dp.subChart !== C) continue;
                    let p=dp.point;
                    p.hide();
					// p.attr({cx:-100,cy:-100});
					moveChartPointTo(p,-100,-100,pointShape);
					if (dp.pointList) {
						for (let j=0; j<dp.pointList.length; j++) {
							dp.pointList[j].hide();
							// dp.pointList[j].attr({cx:-100,cy:-100});
							moveChartPointTo(dp.pointList[j],-100,-100,pointShape);
						}
					}
                    if (p.data('labelObj')) {
                        if (C.dragContext=='pan') {
                            //var tSet=p.data('labelObj');
                            //tSet.hide();
                            //tSet.translate(-p.attr('cx')-50,-p.attr('cy')-50);
                        }
                        else {
                            p.data('labelObj').remove();
                            p.removeData('labelObj');
                        }
                    }
                }
                // curve fitting
                if (mainC.datasets[dsIdx].dataConnection) mainC.datasets[dsIdx].dataConnection.svg.remove();
            }
        }
        if (C.clearChart) {
            C.clearChart();
        }
    };

    var getChartScaleRange = function(optimize,minV,maxV,dimSize) { // dimSize is optional  // *public*
        let deltaV=maxV-minV;
		if (deltaV===0) { // minV=maxV => 1 datapoint
			if (minV) {
				minV-=0.05*minV;
				maxV+=0.05*maxV;
			}
			else { // minV=maxV=0
				minV=-0.05;
				maxV=0.05;
			}
			deltaV=maxV-minV;
        }
        else if (optimize) {
            minV-=(0.05*deltaV); // - 5% range
            maxV=(maxV*1)+(0.05*deltaV); // + 5% range (force to number)
            deltaV*=1.1;
        }
//console.log('Min='+minV+', Max='+maxV+', Delta='+deltaV);
        //Move delta within ]1-10] range
        let shiftFactor=1;
        if (deltaV <= 1) {
            //while (deltaV*shiftFactor < 1) {shiftFactor*=10;}
            let v=Math.floor(1/deltaV)+'';
//console.log(1/deltaV+' => '+v);
            shiftFactor=Math.pow(10,v.length);
        }
        else if (deltaV > 10) {
            //while (deltaV*shiftFactor > 10) {shiftFactor/=10;}
            let v=Math.floor(deltaV - 1)+''; // 100 => 99 (length=2 !3)
//console.log(deltaV+' => '+v);
            shiftFactor=Math.pow(10,-(v.length - 1));
        }
        let optMinV=Math.floor(minV*shiftFactor)/shiftFactor,
        //var maxV=maxV+(0.05*deltaV); // + 5% range
			scaledRange=(maxV-optMinV)*shiftFactor, //deltaV*shiftFactor;
			scaledTick;
        if (dimSize && dimSize < 250) { // 5 ticks
            scaledTick=(scaledRange > 5)? 2 : (scaledRange > 2.5)? 1 : (scaledRange > 1.1)? 0.5 : 0.2;
        }
        else { // 10 ticks
            scaledTick=(scaledRange > 5)? 1 : (scaledRange > 2.5)? 0.5 : (scaledRange > 1.1)? 0.2 : 0.1;
        }
        let tickSize=scaledTick/shiftFactor;
    //console.log('D='+deltaV+', F='+shiftFactor+', minV='+minV+', optMinV='+optMinV+', tick='+tickSize);
        return [minV,optMinV,maxV,tickSize];
    };

	/******************* Point selection & highlight ***********************/
	var setPointSelection = function(p,action,type) { // point click event // *public*
		/* type=max/min. amount ofinfo displayed in popup (only if action=on) */
		let mainC=(p.data('js') && p.data('js').dataset && p.data('js').dataset.params)? p.data('js').dataset.params.chart : null;
//if (mainC && p.data('js').dataset && p.data('js').dataset.params.visible==1) return; // no point selection for dimmed dataset
		if (mainC && mainC.selectAllDatasets) {
			let linkedPoints=mainC.datasetLinks[p.data('js').externalID],
				refAction=(action=='auto' && p.data('showLabel'))? 'off' : (action=='auto')? 'on' : action; // set action based on reference point
			for (let i=0; i < linkedPoints.length; i++) {
				_selectPoint(linkedPoints[i].point,refAction,type);
			}
		}
		else {_selectPoint(p,action,type);}
	};

	function _selectPoint(p,action,type,pX,pY) { // private, optional: type,pX,pY
		if (action=='auto') { // invert selection
			if (p.data('showLabel')) { // already showing label => set OFF
				p.data('showLabel',null);
				_emphasizePoint(p,'off');
			}
			else { // set ON
				p.data('showLabel',type).toFront();
				_emphasizePoint(p,'on',type,pX,pY);
			}
		}
		else if (action=='on') { // set ON
			if (!p.data('showLabel')) {
				p.data('showLabel',type).toFront();
				_emphasizePoint(p,'on',type,pX,pY);
			}
		}
		else { // off
			p.data('showLabel',null);
			_emphasizePoint(p,'off');
		}
		//alert('point ['+p.id+'] '+p.data('index')+' ('+p.data('x')+','+p.data('y')+')');
	}

	var setPointExclusion = function(p) { // *public*
		_selectPoint(p,'off');
		let dp=p.data('js');
		dp.noLine=(dp.noLine)? false : true;
		let dSet=dp.dataset;
		dSet.dataConnection.svg.remove();
		//var pathStrg='',startLine=true;
		let pathPoints=[],
			dpIdx;
		for (let i=0; i < dSet.data.length; i++) {
			let dpi=dSet.data[i],
				pi=dpi.point;
			if (!dpi.noLine) {
				pathPoints.push(pi);
				/*
				//if (startLine==true) {
				//	pathStrg='M'+pi.attr('cx')+','+pi.attr('cy')+' R'; // curve to
				//	startLine=false;
				//}
				//else {pathStrg+=' '+pi.attr('cx')+','+pi.attr('cy');}
				*/
			}
			if (dpi===dp) {dpIdx=i;}
		}
		let mainC=dSet.params.chart;
		if (mainC.connectPoints || dSet.dataConnection) {
			//dSet.line=C.canvas.path(pathStrg).attr({stroke:dSet.params.color,'stroke-width':1}).toBack();
			let C=dp.chart || mainC;
			dSet.dataConnection.svg=drawPath(C,dSet,pathPoints).toBack();
			C.plotArea.toBack();
			C.backPanel.toBack();
		}
	//console.log(dSet.params.index+','+dpIdx);
		if (mainC.onPointExclusion) {mainC.onPointExclusion(dSet.params.index,dpIdx,dp.noLine);}
	};

	var setLabelDisplay = function(p,action,type,e) { // point hover event  // *public*
		let C=p.data('js').dataset.params.chart;
		if (C.datasets.length > 1 && C.datasetLinks && p.data('js').externalID && C.datasetLinks[p.data('js').externalID]) { // points with same id in other datasets, not defined for all charts
			var linkedPoints=C.datasetLinks[p.data('js').externalID];
			for (let i=0; i < linkedPoints.length; i++) {
				let usedType=(C.selectAllDatasets || linkedPoints[i]===p.data('js'))? type : 'none';
				_emphasizePoint(linkedPoints[i].point,action,usedType,null,null,e);
			}
		}
		else {_emphasizePoint(p,action,type,null,null,e);}
	};

	function _emphasizePoint(p,action,type,pX,pY,e) { // private, point hover, click events & search matched
		let pColor,
			dp=p.data('js'),
			mainC=dp.dataset.params.chart,
			C=(dp.subChart)? dp.subChart : mainC;
		if (action==='on') {
			pColor=(mainC.highlightColor)? mainC.highlightColor : idvis.highlightColor;
			//if (p.data('labelObj')) return;
			/* Checking if point is visible (hide if not) */
			if (type !== 'none') {
				let minX=C.plotArea.attr('x'),
					maxX=minX+C.plotArea.attr('width')-1,
					minY=C.plotArea.attr('y'),
					maxY=minY+C.plotArea.attr('height')-1;
				if (!pX && dp.dataset.params.visible) {[pX,pY]=getChartPointCenter(p,dp.dataset.params.pointShape);}
				// if (dp.dataset.params.visible && p.attr('cx') >= minX && p.attr('cx') <= maxX && p.attr('cy') >= minY && p.attr('cy') <= maxY) { // }
				if (dp.dataset.params.visible && pX >= minX && pX <= maxX && pY >= minY && pY <= maxY) {
					if (p.data('labelObj')) {p.data('labelObj').show();}
					else {_displayPointLabel(p,type,e);}
				}
				if (dp.pointList) {
					for (let j=0; j<dp.pointList.length; j++) {
						dp.pointList[j].show();
					}
				}
			}
		}
		else { // off
			if (p.data('showLabel')) { // needed when hovering out of point
				pColor=(dp.dataset.params.visible==1)? idvis.dimColor : (mainC.highlightColor)? mainC.highlightColor : idvis.highlightColor;
			}
			else {
				//pColor=(dp.highlightNames.length)? mainC.highlightedPoints[dp.highlightNames[dp.highlightNames.length-1]].color : dp.dataset.params.color;
				pColor=_findPointColor(dp);
				/*
				pColor=(dp.dataset.params.visible==1)? idvis.dimColor : (dp.color)? dp.color : dp.dataset.params.color;
				if (dp.highlightNames) {
					for (let h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
						if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
							pColor=mainC.highlightedPoints[dp.highlightNames[h]].color;
							break;
						}
					}
				}
				*/
				if (type !== 'none') {
					if (p.data('labelObj')) {
//p.data('labelObj')[1].remove();
						p.data('labelObj').remove();
						p.removeData('labelObj');
					}
					if (dp.pointList) {
						for (let j=0; j<dp.pointList.length; j++) {
							dp.pointList[j].hide();
						}
					}
				}
			}
		}
		p.attr('fill',pColor);
	}

	function _findPointColor(dp,forLabel) { // Checks if point has an highlighting
		let C=(dp.dataset)? dp.dataset.params.chart : null,
			//pColor=(dp.color)? dp.color : (dp.catElement)? dp.catElement.color : (dp.dataset.params.visible==1 && !forLabel)? idvis.dimColor : dp.dataset.params.color;// default
			pColor=(dp.dataset.params.visible==1 && !forLabel)? idvis.dimColor : (dp.color)? dp.color : (dp.catElement)? dp.catElement.color : dp.dataset.params.color;// default
			if (dp.highlightNames && dp.highlightNames.length) {
			for (let h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
				if (C.highlightedPoints[dp.highlightNames[h]].visible) {
					pColor=C.highlightedPoints[dp.highlightNames[h]].color;
					break;
				}
			}
		}
		return pColor;
	}
	
	function _displayPointLabel(p,type,e) { // private, !!! also in peptidePlot.js !!!
		let dp=p.data('js'),
			C=(dp.dataset)? dp.dataset.params.chart : null,
			text=(C && C.customPointLabel)? C.customPointLabel(dp,type) : (dp.info)? dp.info(type) : dp.label, //user-defined text
		//var lColor=(dp.highlightNames.length)? C.highlightedPoints[dp.highlightNames[dp.highlightNames.length-1]].color : dp.dataset.params.color;
			lColor=_findPointColor(dp,true),
			clickData=(C && C.pointOnClick)? [C.pointOnClick,[dp.dataset.params.externalID || dp.dataset.params.index,dp.label,dp.externalID]] : null,
			[pX,pY,pR]=getChartPointCenter(p,dp.dataset.params.pointShape);
		if (clickData && C.categories) {clickData[1].push(C.categories[dp.catIndex].externalID);} // category plot => add category
		//let tSet=drawLabel(p.paper,p,p.attr('cx'),p.attr('cy'),p.attr('r'),p.attr('r'),text,lColor,clickData);
		let tSet=drawLabel(p.paper,p,pX,pY,pR,pR,text,lColor,clickData,e);
		// tSet[0].data({x:p.attr('cx'),y:p.attr('cy')});
		tSet[0].data({x:pX,y:pY});
		p.data({labelObj:tSet});
	}

	var drawLabel = function(paper,refSVG,x,y,d_x,d_y,text,lColor,clickData,e) { // !!! also in peptidePlot.js !!! // *public*
//if (e) console.log('T',e.layerX,e.layerY);
		/* Returns a set containing another set a link line:
		* [[popup_frame,text,sensor],link_to_point]
		*/
		let shift=15,
			t=paper.text(x+d_x+shift,y-d_y-shift,text).attr({'font-size':10,'text-anchor':'start',fill:lColor}), //,'font-weight':'bold','fill-opacity':0.6
			tBox=t.getBBox(),
			tx=tBox.x,
			ty=tBox.y,
			tw=tBox.width,
			th=tBox.height,
		//var tb; //Bubble around text
		//var dy=Math.round((th-6)/2);
			th05=Math.round(th/2);
		if (tx+tw > paper.width) {
			//t.attr({x:paper.width-tBox.width});
			tx=tx-2*(d_x+shift)-tw; //tBox.width;
			t.attr({x:tx});
		}
		if (ty < 0) {
			t.attr({y:2+th05});
			ty=t.getBBox().y;
		}
		if (e && e.layerX >= tx-2 && e.layerX <= tx+tw+2 && e.layerY >= ty-2 && e.layerY <= ty+th+2) { // make sure popup does not cover mouse to prevent immediate hover('off')
//console.log('M',e.layerX,e.layerY);
			t.attr({y:e.layerY+10+th05}); // move 10px down mouse position
			ty=t.getBBox().y;
		}
		let tb=paper.rect(Math.round(tx)-2.5,Math.round(ty)-2.5,tw+5,th+5,4).attr({stroke:lColor,fill:'#fff','fill-opacity':0.9}),
			tl=_drawLabelLink(paper,x,y,d_x,d_y,tb);
		t.toFront();
		let ts=paper.rect(Math.round(tx)-2.5,Math.round(ty)-2.5,tw+5,th+5,4).attr({stroke:'none',fill:'#fff','fill-opacity':0}), // sensor
			tSet=paper.set(paper.set(tb,t,ts),tl),
			popupCoord=null,
			cursorType=(clickData)? 'pointer' : 'default';
		ts.attr({cursor:cursorType}).data({set:tSet,dragContext:false}) // js:p,
		.hover(
			function() { // in
				let tSet=this.data('set');
				if (!tSet) return; // hovering is faster that tSet creation/registration
				tSet[0][0].attr({stroke:idvis.highlightColor}); tSet[1].attr({stroke:idvis.highlightColor}); // frame & connecting line
				if (idvis.closePopupIcon) {idvis.closePopupIcon.remove(); idvis.closePopupIcon=null;} // already drawn for another popup
				if (this.data('dragContext')===true) return;
				let clx=this.attr('x')+this.attr('width')-13, cly=this.attr('y')+3;
				idvis.closePopupIcon=paper.set(paper.rect(clx,cly,10,10,2).attr({stroke:'none',fill:'#f00'}),
											   paper.path('M,'+(clx+2)+','+(cly+2)+' l6,6 M'+(clx+2)+','+(cly+8)+'l6,-6').attr({stroke:'#fff','stroke-width':2})
				).attr({cursor:'pointer'}).data({'set':tSet})
				.hover(
					function() {
						let tSet=this.data('set');
						tSet[0][0].attr({stroke:idvis.highlightColor}); tSet[1].attr({stroke:idvis.highlightColor}); // frame & connecting line
						idvis.closePopupIcon.show();
					},
					function() {}
				)
				.click(
					function() {
						idvis.closePopupIcon.remove();
						idvis.closePopupIcon=null;
						let closeCallback=(refSVG.data('js').select)? refSVG.data('js').select : setPointSelection; 
						closeCallback(refSVG,'off');
					} // TODO: Check popup for other elements than points!!!!!
				);
			},
			function() { // out
				let tSet=this.data('set');
				if (!tSet) return; // hovering is faster that tSet creation/registration
				let color=tSet[0][1].attr('fill'); // original color
				tSet[0][0].attr({stroke:color}); tSet[1].attr({stroke:color}); // frame & connecting line
				if (idvis.closePopupIcon) {idvis.closePopupIcon.hide();} // hide, not remove in case hover moved to closePopupIcon itself!!
			}
		)
		.drag( // this = tset[0]
			function(dx,dy) { // extend
				let tSet=this.data('set');
				if (this.data('dragContext')===false) { // only for 1st dx/dy
					this.data({dragContext:true});
					tSet[1].remove(); // remove connecting line
					if (idvis.closePopupIcon) idvis.closePopupIcon.remove();
					idvis.closePopupIcon=null;
				}
				tSet[0][0].attr({x:popupCoord[0]+dx,y:popupCoord[1]+dy}); // frame
				tSet[0][1].attr({x:popupCoord[2]+dx,y:popupCoord[3]+dy}); // text
				tSet[0][2].attr({x:popupCoord[0]+dx,y:popupCoord[1]+dy}); // sensor
			},
			function() { // start (also triggered by simple click)
				let tSet=this.data('set');
				this.attr({cursor:'move'});
				popupCoord=[tSet[0][0].attr('x'),tSet[0][0].attr('y'),tSet[0][1].attr('x'),tSet[0][1].attr('y')];
				//p.data('labelObj')[1].remove(); // remove connecting line
			},
			function() { // end (also triggered by simple click)
				if (this.data('dragContext')===true) { // not a simple click
					if (!clickData) this.data({dragContext:false});
					let tSet=this.data('set');
					this.attr({cursor:cursorType});
					/* redraw & store connecting line */
					let tSet2=paper.set(tSet[0],_drawLabelLink(paper,x,y,d_x,d_y,tSet[0][0]));
					this.data({set:tSet2}); // update to new set
//console.log('DRAG',this.data());
					//refSVG.removeData('labelObj'); // Clear ... necessary?
					refSVG.data({labelObj:tSet2}); // & rebuild is necessary
				}
			}
		);
		if (clickData) {
//console.log('OK');
			ts.click(function(){
				if (this.data('dragContext')===true) {this.data({dragContext:false});}
				//else {clickData[0](clickData[1],clickData[2],clickData[3]);}  // dp.dataset.params.index,dp.label,dp.externalID
				else {clickData[0](...clickData[1]);}  // <dataset id|index>,<element label>,<element externalID>[,<category id|index>]
			});
		}

		return tSet;
	};

	function _drawLabelLink(paper,x,y,d_x,d_y,tb) {
		//var shift=15;
		let tx=tb.attr('x'),
			ty=tb.attr('y'),
			tw=tb.attr('width'),
			th=tb.attr('height'),
			th05=Math.round(th/2),
			txw=tx+tw,
			tyh=ty+th,
			sx;
//console.log(x,y,document.elementsFromPoint(x,y));
//d_x=d_y=0;
		if (Math.abs(tx-x) <= Math.abs(txw-x)) {sx=tx+th05;} //point closer to popup left side
		else {sx=txw-th05;} //point closer to popup right side

		let sy=ty+th05,
			tl=paper.path('M'+sx+' '+sy+'L'+x+' '+y).attr({stroke:tb.attr('stroke')}); // popup->point // ,'stroke-width':2,'stroke-opacity':0.4
		if (Math.abs(sx-x) >= Math.abs(sy-y)) { // dx >= dy
			if (x<tx) { // popup on the right
//console.log(1,1);
				tl.attr({'clip-rect':(x+d_x)+',0,'+(tx-x-d_x)+','+paper.height});
//paper.rect(x+d_y,0,tx-x-d_x,paper.height,0);
			}
			else if (x>txw) {  // popup on the left
//console.log(1,2);
				tl.attr({'clip-rect':txw+',0,'+(x-d_x-txw)+','+paper.height});
//paper.rect(txw,0,x-d_x-txw,paper.height,0);
			}
		}
		else { // dy > dx
			if (y<ty) { // popup below
//console.log(2,1);
				tl.attr({'clip-rect':'0,'+(y+d_y)+','+paper.width+','+(ty-y-d_y)});
//paper.rect(0,y+d_y,paper.width,ty-y-d_y,0);
			}
			else if (y>tyh) { // popup above
//console.log(2,2);
				tl.attr({'clip-rect':'0,'+tyh+','+paper.width+','+(y-d_y-tyh)});
//paper.rect(0,tyh,paper.width,y-d_y-tyh,0);
			}
		}
//console.log(tl.attr('clip-rect'));

		return tl;
	}

	function _changeLabelColor(labelObj,newColor) { // labelObj is paper.set(paper.set(tb,t,ts),tl),
		labelObj[0][0].attr('stroke',newColor); // box
		labelObj[0][1].attr('fill',newColor); // text
		labelObj[1].attr('stroke',newColor); // link
	}
	
	function _moveLabelsToFront(mainC) { // Moves all displayed labels to front to make sure they stay visible after a search or highlighting
		for (let dsIdx=0; dsIdx<mainC.datasets.length; dsIdx++) {
			for (let i=0; i < mainC.datasets[dsIdx].data.length; i++) {
				let p=mainC.datasets[dsIdx].data[i].point;
				if (p.data('labelObj')) {
					p.data('labelObj').toFront();
				}
			}
		}
		_moveDragAreasToFront(mainC);
	}
	function _moveDragAreasToFront(mainC) {
		if (mainC.subChart) { // multi-chart
			for (let c=0; c < mainC.activeCharts.length; c++) {
				if (mainC.subChart[mainC.activeCharts[c]].dragArea) {mainC.subChart[mainC.activeCharts[c]].dragArea.toFront();}
			}
		}
		else if (mainC.dragArea) {mainC.dragArea.toFront();}
	}

	/********** Path drawing **********/
	var drawPath = function(C,element,pathPoints) { // public
//console.log((element.params)? element : 'ok');
		let [connectObj,pointShape]=(element.params)? [element.dataConnection,element.params.pointShape] : [element,'circle'] // element: dataset or feature
		    type=connectObj.lineType || 'line', // default
			color=connectObj.color || '#000',
			pattern=connectObj.pattern || null,
			width=connectObj.width || 1,
			opacity=connectObj.opacity || 1;
		if (element.params && element.params.visible==1) { // dimmed dataset
			color=idvis.dimColor;
			opacity=idvis.dimSetOpacity;
		}
/*
		if (params.isFeature) { // feature (threshold)
			type=params.lineType;
			color=params.color;
			pattern=params.pattern;
			width=1;
			opacity=params.opacity;
		}
		else if (params.line) {
			if (params.line.type) type=params.line.type;
			if (params.line.color) color=params.line.color;
			if (params.line.pattern) pattern=params.line.pattern;
			if (params.line.width) width=params.line.width;
			if (params.line.opacity) opacity=params.line.opacity;
		}
*/		
		
		let pathArray=[];
		if (!type.match('line|regression') && pathPoints.length <= 2) type='line';
		let prevGap=true,
		    numConnPoints=0,
			lCode=(type=='curve')? 'R' : 'L';
		let i;
		for (i=0; i<pathPoints.length; i++) {
			if (!pathPoints[i]) { // Gap in path
				if (!type.match('line|regression') && numConnPoints == 2) { // convert prev section to line if <= 2 points
					let j=pathArray.length-1-numConnPoints;
					if (j > 0 && pathArray[j] != 'M') pathArray[j]='L';
//console.log(j,numConnPoints,pathArray);
				}
				prevGap=true;
				numConnPoints=0;
				if (!idvis.isNumber(pathArray[pathArray.length-1])) pathArray.pop();
				continue;
			}
			let x,y;
			if (connectObj.isFeature || type=='regression') {[x,y]=pathPoints[i];} //{x=pathPoints[i][0]; y=pathPoints[i][1];}
			// else {x=pathPoints[i].attr('cx'); y=pathPoints[i].attr('cy');}
			else {[x,y]=getChartPointCenter(pathPoints[i],pointShape);}
			if (prevGap) { //if (i===0) {
				pathArray.push('M',x,y);
				if (i < pathPoints.length-1) pathArray.push(lCode);
			}
			else {pathArray.push(x,y);}
			prevGap=false;
			numConnPoints++;
		}
		if (!type.match('line|regression') && numConnPoints == 2) { // convert prev section to line if <= 2 points
			let j=pathArray.length-1-numConnPoints;
			if (j > 0 && pathArray[j] != 'M') pathArray[j]='L';
//console.log(j,numConnPoints,pathArray);
		}		
		let pathStrg=pathArray.join(' ');
//console.log(pathStrg);
		let clipStrg=C.plotArea.attr('x')+','+C.plotArea.attr('y')+','+C.plotArea.attr('width')+','+C.plotArea.attr('height'),
		    paper=(C.mainChart)? C.mainChart.canvas : C.canvas,
			path=paper.path(pathStrg).attr({stroke:color,'stroke-width':width,'stroke-opacity':opacity,'clip-rect':clipStrg});
		//??? for (let i=0; i<pathPoints.length; i++) {pathPoints[i].toFront();}
		if (pattern) path.attr({'stroke-dasharray':pattern});
		return path;
	};

	/******************* Mouse drag management *********************/
	var setDragging = function(C,x,y,e) { // x,y are absolute to document // *public*
		//var mouseInChart=getMousePositionInChart(e);
		let mousePos=idvis.getMousePositionInElement(e),
//debug('ctrl='+idvis.modKeyPressed+' : '+x+'('+mousePos[0]+'), '+y+'('+mousePos[1]+')');
			mouseX=mousePos[0]+C.plotArea.attr('x')-1,mouseY=mousePos[1]+C.plotArea.attr('y')-1;
		if (idvis.modKeyPressed) {
			C.dragContext='pan';
			C.plotArea.attr({cursor:'move'});
			C.panProcess.x=0;//mouseX;
			C.panProcess.y=0;//mouseY;
			C.panProcess.dx=0;
			C.panProcess.dy=0;
		}
		else {
			C.dragContext='area';
			C.plotArea.attr({cursor:'crosshair'});
			//_setDragArea(C,mouseInChart[0],mouseInChart[1]);
			_setDragArea(C,mouseX,mouseY);
		}
	};

	var extendDragging = function(C,dx,dy) { // ...,x,y,e *public*
		if (C.dragContext=='pan') {
			C.panProcess.dx=dx-C.panProcess.x;
			C.panProcess.dy=dy-C.panProcess.y;
//console.log('Pan: x='+VP.panProcess.x+', y='+VP.panProcess.y+', dx='+VP.panProcess.dx+', dy='+VP.panProcess.dy);
			if (Math.abs(C.panProcess.dx) > 5 || Math.abs(C.panProcess.dy) > 5) {
				panChart(C);
				C.panProcess.x+=C.panProcess.dx;
				C.panProcess.y+=C.panProcess.dy;
			}
		}
		else {_extendDragArea(C,dx,dy);}
	};

	var endDragging = function(C) { // *public*
		if (C.dragContext=='pan') {
			C.dragContext='area';
			C.panProcess.x=0;
			C.panProcess.y=0;
			//panChart(C);
			C.panProcess.dx=0;
			C.panProcess.dy=0;
//console.log('END');
			/*** Clear & plot chart ***/
			clearChart(C);
			plotChart(C);
		}
		else {
			_endDragArea(C);
			//C.dragContext=null;
		}
		C.plotArea.attr({cursor:'auto'});
	};

	/******** Drag area management *******/
	function _setDragArea(C,mX,mY) { // private, position in chart
		let mainC=(C.mainChart)? C.mainChart : C,
			x,y,w,h;
		if (C.noZoomX && !mainC.datasets) {x=C.plotArea.attr('x'); w=C.plotArea.attr('width');} //
		else {x=mX; w=0;}
		if (C.noZoomY && !mainC.datasets) {y=C.plotArea.attr('y'); h=C.plotArea.attr('height');} //
		else {y=mY; h=0;}
		if (C.dragArea.data('status') == 'on') {
			C.dragArea.attr({width:w,height:h});
		}
		C.dragArea.data({startX:x,startY:y}) //,status:'extend'
		.attr({x:x,y:y,width:w,height:h})
		.toFront()
		.show();
//console.log('SET:'+x+','+y+','+w+','+h);
	}

	function _extendDragArea(C,dx,dy) { // private
		let mainC=(C.mainChart)? C.mainChart : C,
			dragArea=C.dragArea,
			plotArea=C.plotArea,
			chartX0=plotArea.attr('x')-1,
			chartY0=plotArea.attr('y')-1;
			if (!C.noZoomX || mainC.datasets) {
			if (dx > 0) {
				dragArea.attr({width:dx});
				if (dragArea.attr('x')+dragArea.attr('width') > plotArea.attr('width')+chartX0) {dragArea.attr({width:plotArea.attr('width')-dragArea.attr('x')+chartX0});}
			}
			else {
				dragArea.attr({x:dragArea.data('startX')+dx,width:-dx});
				if (dragArea.attr('x') < plotArea.attr('x')) {
					let extra=plotArea.attr('x') - dragArea.attr('x');
					dragArea.attr({x:plotArea.attr('x'),width:dragArea.attr('width')-extra});
				}
			}
		}
		if (!C.noZoomY || mainC.datasets) {
			if (dy > 0) {
				dragArea.attr({height:dy});
				if (dragArea.attr('y')+dragArea.attr('height') > plotArea.attr('height')+chartY0) {dragArea.attr({height:plotArea.attr('height')-dragArea.attr('y')+chartY0});}
			}
			else {
				dragArea.attr({y:dragArea.data('startY')+dy,height:-dy});
				if (dragArea.attr('y') < plotArea.attr('y')) {
					let extra=plotArea.attr('y') - dragArea.attr('y');
					dragArea.attr({y:plotArea.attr('y'),height:dragArea.attr('height')-extra});
				}
			}
		}
//console.log('EXT:'+dragArea.attr('x')+','+dragArea.attr('y')+','+dragArea.attr('width')+','+dragArea.attr('height'));
	}

	function _endDragArea(C) { // private
		let mainC=(C.mainChart)? C.mainChart : C;
		if ((!mainC.datasets && ((C.noZoomX && C.dragArea.attr('height') > 5) || (C.noZoomY && C.dragArea.attr('width') > 5))) || ((mainC.datasets || !C.noZoomX && !C.noZoomY) && C.dragArea.attr('width')*C.dragArea.attr('height') > 25)) { // big enough => end extension
			C.dragArea.data({status:'on'});
			if (C.autoZoom) {
				zoomIn(C);
			}
		}
		else { // too small => hide
			_clearDragArea(C.dragArea);
		}
	}
	
	function _clearDragArea(dragArea) {
		dragArea.hide()
		.data({status:'off'})
		.attr({width:0,height:0});
	}

	/************* Chart Export *************/
	var exportSVGtoImg = function(svgDivId,imgName,exportScript,format='png') {
		//if (!format) format='png';

		//put html svg to svgHTML
		//svgFix make svg standard
		let svgHTML = document.getElementById(svgDivId).innerHTML,
			data;
		if (format==='png') { // specific processing for png using canvas element
			//create canvas element
			let canvas = document.createElement('canvas');
			canvas.id = svgDivId+"_exportImg";
			document.body.appendChild(canvas);
			//canvas.style.display='none';

			//canvg is a SVG parser and renderer.
			//It takes a URL to a SVG file or the text of an SVG file,
			//parses it in JavaScript, and renders the result on a Canvas element
			canvg(canvas, svgfix(svgHTML),{
				renderCallback:function() { // Called when image load is complete!!!
						//Returns the content of the current canvas as an image (data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAATYAAA......)
						data = canvas.toDataURL("image/png");
//console.log(data);
						//transform data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAATYAAA.... to iVBORw0KGgoAAAANSUhEUgAAATYAAA....
						data = data.substr(data.indexOf(',') + 1).toString();
						document.body.removeChild(canvas);
						_exportDataStringImage();
					}
				}
			);
		}
		else if (format=='svg') {
			data = svgfix(svgHTML);
			_exportDataStringImage();
		}

		function _exportDataStringImage() { // private
			//create temporary input type hidden with image
			var dataInput = document.createElement("input");
			dataInput.setAttribute("name", 'imgdata');
			dataInput.setAttribute("value", data);
			dataInput.setAttribute("type", "hidden");

			//create temporary input type with name of image
			var nameInput = document.createElement("input") ;
			nameInput.setAttribute("name", 'name');
			nameInput.setAttribute("value", imgName+'.'+format);
			nameInput.setAttribute("type", "hidden");

			//create temporary html form, update body, submit and remove html form
			var myForm = document.createElement("form");
			myForm.method = 'post';
			myForm.action = exportScript; //"./exportSVG.cgi";
			myForm.appendChild(dataInput);
			myForm.appendChild(nameInput);
			document.body.appendChild(myForm);
			myForm.submit();
			document.body.removeChild(myForm);
		}
	};
	
	/************* Linear Regression *************/
	// https://stackoverflow.com/questions/6195335/linear-regression-in-javascript
	var linearRegression = function(x,y) {
        var lr = {},
            n = y.length,
            sum_x = 0,
            sum_y = 0,
            sum_xy = 0,
            sum_xx = 0,
            sum_yy = 0;
		if (n==0) return lr;
		
        for (let i = 0; i < n; i++) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += (x[i]*y[i]);
            sum_xx += (x[i]*x[i]);
            sum_yy += (y[i]*y[i]);
        } 

        lr.slope = (n*sum_xy - sum_x*sum_y) / (n*sum_xx - sum_x*sum_x);
        lr.intercept = (sum_y - lr.slope*sum_x)/n;
        lr.r = (n*sum_xy - sum_x*sum_y)/Math.sqrt((n*sum_xx-sum_x*sum_x)*(n*sum_yy-sum_y*sum_y));
		lr.r2 = Math.pow(lr.r,2);

        return lr;
	};

	/****** return idvis public methods *****/
	return {
		registerChart : registerChart,
		addDataset : addDataset, addDataSet : addDataset, // just to be safe
		generateDatasetAES : generateDatasetAES,
		drawChartLayout : drawChartLayout,
		initializeForm : initializeForm,
		initializeChart : initializeChart,
		computeDefaultRange : computeDefaultRange,
		plotChart : plotChart,
		drawChartPoint : drawChartPoint,
		drawFeature : drawFeature,
		drawPath : drawPath,
		setDatasetDisplay : setDatasetDisplay,
		addHighlighting : addHighlighting,
		deleteHighlighting : deleteHighlighting,
        editHighlighting : editHighlighting,
        applyHighlightEdition : applyHighlightEdition,
        selectHighlightedPoints : selectHighlightedPoints,
        setHighlightVisibility : setHighlightVisibility,
        setUnmatchedPointsVisibility : setUnmatchedPointsVisibility,
        searchDataPoints : searchDataPoints,
        manageSelection : manageSelection,
        listSelectedPoints : listSelectedPoints,
        selectThreshold : selectThreshold,
        updateThreshold : updateThreshold,
		activateDimArea : activateDimArea,
		toggleColorSelection : toggleColorSelection,
		changeElementColor : changeElementColor,
		zoomIn : zoomIn,
		zoomOut : zoomOut,
        panChart : panChart,
        clearChart : clearChart,
        getChartScaleRange : getChartScaleRange,
		setPointSelection : setPointSelection,
		setPointExclusion : setPointExclusion,
		setLabelDisplay : setLabelDisplay,
		drawLabel : drawLabel,
		setDragging : setDragging,
		extendDragging : extendDragging,
		endDragging : endDragging,
		exportSVGtoImg : exportSVGtoImg,
		linearRegression : linearRegression
	};
})();


/*
####>Revision history<####
# 2.2.0 [UPDATE] Point shapes support, full color editing, form positioning & multiple minor improuvements (PP 06/05/21) 
# 2.1.1 [FEATURE] Partial support for interactive color management of datasets (PP 03/12/20)
# 2.1.0 [FEATURE] Linear regression & custom elements in main form  & improved highliting display (PP 01/09/20)
# 2.0.4 [FEATURE] Handles an optional external search-text conversion function (PP 09/12/19)
# 2.0.3 Double click removes drag area on non-zoomable charts (PP 30/05/19)
# 2.0.2 [Fix] bug in feature function min/max Y computation & in clear when multi-value data points & [Update] SVG image export (PP 06/06/18)
# 2.0.1 [Update] Checks for callback 'onThresholdUpdate' function after threhold edition & [fix] label popup link after drag (PP 29/09/17)
# 2.0.0 Updated chartLibrary2.js to code-encapsulated idvis.js (PP 13/06/17)
# 1.4.3 Improved chart.showDatasets detection & chart-specific form elements (PP 20/11/16)
# 1.4.2 Bug fix in new pointOnList behavior (PP 10/11/16)
# 1.4.1 pointOnList now also accepts an array of [action,function] array (PP 03/11/16)
# 1.4.0 Drawing legends for highlights (PP 23/10/16)
# 1.3.8 Better handling of single data poin in getChartScaleRange (PP 27/05/16)
# 1.3.7 Added chart registering (PP 20/02/16)
# 1.3.6 max-height on Datasets FIELDSET (PP 28/10/15)
# 1.3.5 Handle dataset.params.pointSizeRule for scatter plots (PP 21/10/15)
# 1.3.4 Increase maximum search matches to 100 (PP 12/10/15)
# 1.3.3 New getMousePositionInElement() to fix bug in mouse position relative to parent element (PP 24/07/15)
# 1.3.2 Conditional 'stroke-dasharray' attribute on threshold lines to prevent invisibility on some systems (PP 22/04/15)
# 1.3.1 ThresholdLine edition now calls moveThresholdLine() instead of plotChart() to redraw line (PP 21/04/15)
# 1.3.0 Upgrade addHighlighting function to allow partial ID match. Uses pattern eg. '^###-' as optional 4th parameter (PP 07/11/14)
# 1.2.9 Improved mouse position detection for IE 9+ & canvas export (PP 09/10/14)
# 1.2.8 getMousePositionInChart() non longer used (PP 11/08/14)
# 1.2.7 Minor behavior bug fix in Hide/show points not matching highlight & string shortening function (PP 08/07/14)
# 1.2.6 Handles svg export as svg (PP 24/06/14)
# 1.2.5 Hide/show points not matching highlight (PP 15/06/14)
# 1.2.4 Compatibility with spectrumPlot.js (PP 01/04/14)
# 1.2.3 Added exportSVGtoImg function & conditional export button (PP 11/02/14)
# 1.2.2 GPL license (PP 23/09/13)
# 1.2.1 Added point opacity management (PP 15/07/13)
# 1.2.0 Handles sub graphes (PP 17/03/13)
# 1.1.0 Better check on point linkage across datasets during mouseover (PP 29/01/13)
# 1.0.9 Added startValue to threshold object (PP 02/01/13)
# 1.0.8 Tests if point exclusion if allowed (PP 05/12/12)
# 1.0.7 Uses shift, crtl & alt as modifier keys (PP 03/12/12)
# 1.0.6 Compatibility with point exclusion from line plot (PP 02/12/12)
# 1.0.5 Form position option (PP 30/05/12)
# 1.0.4 Fix bug single dataPoint (PP 24/05/12)
# 1.0.3 Minor changes (PP 13/05/12)
# 1.0.2 Compatibility with line plot (PP 12/05/12)
# 1.0.1 Check value provided for threshold change (PP 04/05/12)
# 1.0.0 Production version
# 0.9.9 Fix bug number of points affected when renaming highlighting (PP 29/03/12)
*/
