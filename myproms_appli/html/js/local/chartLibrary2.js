/*
################################################################################
# chartLibrary2.js         1.5.0                                               #
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

// future name: biovis.js
var chartLibrary=true; // global variable to indicate the library is used
var cubiojsCharts=[]; // Global: needed to call chart functions from its own dynamically generated html
var modKeyPressed=false; // global


/*****************************************************
  Overwrites Raphaël text to fix bug in Chrome when drawing text in hidden (display:none) div => Generates collateral errors in other browsers
******************************************************/
//(function() { //
//	var oldText = Raphael.prototype.text;
//
//	Raphael.prototype.text = function () {
//		var textElement = oldText.apply(this, arguments);
//
//		setTimeout(function () {
//			var tspanElement = textElement.node.getElementsByTagName('tspan')[0];
//			if (tspanElement) tspanElement.setAttribute('dy', 0);
//		}, 1);
//
//		return textElement;
//	};
//})();


/*****************************************************
                   Generic Objects
******************************************************/

/************************ Threshold line object ********************/
function thresholdLine(C,axis,name,initValue,color,dashString,keepAbove1,roundValue,editable) {
	this.chart=C;
	this.axis=axis;
	this.name=name;
	this.color=(color)? color : '#000';
	this.keepAbove1=(keepAbove1)? true :false;
	this.roundValue=roundValue;
	this.path=null;
	this.pathPos=null;
	this.popup=null;
	this.startValue=initValue;
	this.setValues(initValue);
	this.dash=(dashString)? dashString : null;
	this.editable=editable;
}
thresholdLine.prototype={
	setValues: function(initValue) {
		this.initValue=initValue;
		this.value=(this.chart.convertValue)? this.chart.convertValue(this.axis,initValue) : initValue;
		if (this.name) {
			var dispValue;
			if (initValue < 1 && this.keepAbove1) {
				dispValue=(this.roundValue)? '1/'+(1/initValue).toFixed(this.roundValue) : '1/'+(1/initValue);
			}
			else {
				dispValue=(this.roundValue)? (initValue*1).toFixed(this.roundValue) : initValue;
			}
			this.label=this.name+': '+dispValue;
		}
		else {this.label=null;}
	}
}


/*****************************************************
                   Generic Functions
******************************************************/

/******************* Register chart into cubiojs array *********************/
function registerInCuBioJS(Chart) { // private
	cubiojsCharts.push(Chart);
	return cubiojsCharts.length-1; // index of chart in array
}

/******************* HTML interface *********************/
function drawChartBorders(mainC,canvasW,canvasDivID,formDivID) { // *public*
	var borderHTML;
	if (mainC.formPosition=='top') {
		borderHTML='<TABLE><TR><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD></TR></TABLE>';
	}
	else if (mainC.formPosition=='bottom') {
		borderHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD></TR><TR><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
	}
	else if (mainC.formPosition=='left') {
		borderHTML='<TABLE><TR><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD></TR></TABLE>';
	}
	else { // right
		borderHTML='<TABLE><TR><TD style="min-width:'+canvasW+'px"><DIV id="'+canvasDivID+'"></DIV></TD><TD valign="top" style="min-width:200px"><DIV id="'+formDivID+'"></DIV></TD></TR></TABLE>';
	}
	mainC.div.innerHTML=borderHTML;
}
function initializeForm(mainC) { // *public*
	if (mainC.chartID==undefined) { // just to be safe
		mainC.chartID=registerInCuBioJS(mainC);
	}
	var htmlString='';

	/* DataSets visibility */
	if (mainC.showDataSets || (mainC.showDataSets != false && mainC.dataSets.length > 1)) {
		htmlString+='<FIELDSET style="padding:2px;white-space:nowrap;max-height:200px;overflow:auto"><LEGEND><B>Datasets:</B></LEGEND>';
		for (var i=0; i<mainC.dataSets.length; i++) {
			if (i > 0) {htmlString+='<BR>';}
			htmlString+='<INPUT type="checkbox" value="'+i+'" onclick="setDataSetDisplay(cubiojsCharts['+mainC.chartID+'],this.value,this.checked)" checked><FONT style="color:'+mainC.dataSets[i].params.color+'">'+mainC.dataSets[i].params.name+'</FONT>';
		}
		htmlString+="</FIELDSET>\n";
	}

	/* Selection buttons */
	if (mainC.selectable) {
		htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Selection:</B></LEGEND>';
		if (mainC.dataSets.length > 1) {htmlString+='<INPUT type="checkbox" value="1" onclick="cubiojsCharts['+mainC.chartID+'].selectAllDataSets=this.checked"><FONT style="font-weight:bold;font-size:12px">Extend to all datasets</FONT><BR>';}
		htmlString+='<INPUT type="button" value="Select" style="font-weight:bold;font-size:12px;width:90px" onclick="manageSelection(cubiojsCharts['+mainC.chartID+'],\'on\')">';
		htmlString+='<INPUT type="button" value="Un-Select" style="font-weight:bold;font-size:12px;width:90px" onclick="manageSelection(cubiojsCharts['+mainC.chartID+'],\'off\')">';
		if (mainC.pointOnList) {
			if (typeof(mainC.pointOnList)=='function') { // default
				htmlString+='<BR><INPUT type="button" value="List selected" style="font-weight:bold;font-size:12px;width:180px" onclick="listSelectedPoints(cubiojsCharts['+mainC.chartID+'],-1)">';
			}
			else { // array of ['action name',function] arrays
				for (var i=0; i<mainC.pointOnList.length; i++) {
					htmlString+='<BR><INPUT type="button" value="'+mainC.pointOnList[i][0]+' selected" style="font-weight:bold;font-size:12px;width:180px" onclick="listSelectedPoints(cubiojsCharts['+mainC.chartID+'],'+i+')">';
				}
			}
		}
		htmlString+="</FIELDSET>\n";
	}

	/* Search */
	if (mainC.searchable) {
		htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Search:</B></LEGEND>';
		var searchBoxID=mainC.divID+'_search';
		var searchResDivID=mainC.divID+'_srchRes';
		htmlString+='<INPUT type="text" id="'+searchBoxID+'" style="width:140px" value=""><INPUT type="button" value="Go" style="font-weight:bold;font-size:12px;width:40px" onclick="searchDataPoints(cubiojsCharts['+mainC.chartID+'],document.getElementById(\''+searchBoxID+'\').value,\''+searchResDivID+'\',';
		if (typeof(mainC.searchable)=='object') {
			var extSearchText=mainC.searchable.text || 'Extended search',
			    extSearchChkID=mainC.divID+'_extSearch';
			htmlString+='\''+extSearchChkID+'\')"><BR><INPUT type="checkbox" id="'+extSearchChkID+'" value="1"><LABEL for="'+extSearchChkID+'">'+extSearchText+'</LABEL>';
		}
		else {htmlString+='null)">';}
		htmlString+='<DIV id="'+searchResDivID+'"></DIV>'
		htmlString+="</FIELDSET>\n";
	}

	/* Thresholds control */
	if (mainC.editThreshold && mainC.thresholdLines.length) {
		htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Thresholds:</B></LEGEND>';
		var thresSelID=mainC.divID+'_thSel';
		var thresBoxID=mainC.divID+'_thBox';
		if (mainC.thresholdLines.length>1) {
			htmlString+='<SELECT id="'+thresSelID+'" style="font-size:12px;width:180px" onchange="selectThreshold(cubiojsCharts['+mainC.chartID+'],this.value,\''+thresBoxID+'\')"><OPTION value="-1">-= Select =-</OPTION>';
			for (var i=0; i<mainC.thresholdLines.length; i++) {
				htmlString+='<OPTION value="'+i+'">'+mainC.thresholdLines[i].name+'</OPTION>';
			}
			htmlString+='</SELECT>';
		}
		else {htmlString+='<FONT style="font-size:12px;">'+mainC.thresholdLines[0].name+'</FONT>';}
		htmlString+='<BR><FONT style="font-size:12px;">Value: </FONT><INPUT type="text" id="'+thresBoxID+'" style="font-size:12px;width:60px" value=""><INPUT type="button" value="Go" style="font-weight:bold;font-size:12px;width:40px" onclick="updateThreshold(cubiojsCharts['+mainC.chartID+'],\''+thresSelID+'\',\''+thresBoxID+'\')">';
		htmlString+="</FIELDSET>\n";
	}

	/* Point highlighting */
	if (mainC.allowHighlight) {
		var highlightDivID=mainC.divID+'_hlight';
		htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Highlighting:</B></LEGEND><DIV id="'+highlightDivID+'">None</DIV></FIELDSET>\n';
	}

	/* Chart-specific form elements */
	if (mainC.addToForm) htmlString+=mainC.addToForm();

	/* Export image button */
	if (mainC.exportAsImage) {
		if (mainC.exportAsImage[3]) { // image format
			htmlString+='<BR><INPUT type="button" value="'+mainC.exportAsImage[0]+'" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg(\''+mainC.getDivID()+'\',\''+mainC.exportAsImage[1]+'\',\''+mainC.exportAsImage[2]+'\',\''+mainC.exportAsImage[3]+'\')"/>\n';
		}
		else {
			htmlString+='<BR><FONT style="font-weight:bold;font-size:10px">'+mainC.exportAsImage[0]+':</FONT>';
			//png
			htmlString+='<INPUT type="button" value="PNG" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg(\''+mainC.getDivID()+'\',\''+mainC.exportAsImage[1]+'\',\''+mainC.exportAsImage[2]+'\',\'png\')"/>';
			//svg
			htmlString+='<INPUT type="button" value="SVG" style="font-weight:bold;font-size:10px" onclick="exportSVGtoImg(\''+mainC.getDivID()+'\',\''+mainC.exportAsImage[1]+'\',\''+mainC.exportAsImage[2]+'\',\'svg\')"/>\n';

		}
	}

	mainC.formDiv.innerHTML=htmlString;
}

/******************* Chart drawing *********************/
function initializeChart(mainC) { // *public*
	var canvas=mainC.canvas;

	var chartList=[];
	if (mainC.subChart) {
		for (var c=0; c<mainC.activeCharts.length; c++) {
			chartList.push(mainC.subChart[mainC.activeCharts[c]]);
		}
	}
	else {chartList.push(mainC);}

	// Looping through displayed charts
	for (var c=0; c<chartList.length; c++) {
		var C=chartList[c];
//C.plotArea.attr({stroke:'#F00'});
		var chartX=C.plotArea.attr('x')-0.5; // var chartX0=chartX-1;
		var chartY=C.plotArea.attr('y')-1.5;
		var chartW=C.plotArea.attr('width')+2;
		var chartH=C.plotArea.attr('height')+2;

		/* Axes */
		C.usedAxes={};
		var axes=['X','Y','X2','Y2'];
		for (var i=0; i<axes.length; i++) {if (C['axis'+axes[i]+'text'] != null) {C.usedAxes[axes[i]]=true;}}// axis is declared
//console.log(C.usedAxes);
		if (mainC.axisClosure) {canvas.path('M'+chartX+' '+chartY+' l0 '+chartH+' l'+chartW+' 0 l0 -'+chartH+' Z');}
		else {
			if (C.usedAxes.X && !C.noAxisX) {canvas.path('M'+chartX+' '+(chartY+chartH)+' l'+chartW+' 0');} // X
			if (C.usedAxes.Y && !C.noAxisY) {canvas.path('M'+chartX+' '+chartY+' l0 '+chartH);} // Y
			if (C.usedAxes.X2 && !C.noAxisX2) {canvas.path('M'+chartX+' '+chartY+' l'+chartW+' 0');} // X2
			if (C.usedAxes.Y2 && !C.noAxisY2) {canvas.path('M'+(chartX+chartW)+' '+chartY+' l0 '+chartH);} // Y2
		}
		if (C.axisXtext) {
			var titleX=canvas.text(chartX+chartW/2,chartY+chartH+30,C.axisXtext).attr({'font-weight':'bold','font-size':14});
			if (C.axisXtitle) C.axisXtitle=titleX; // used by pcaPlot.js to redraw title
		}
		if (C.axisYtext) {
			var tx=chartX-45; //Math.max(10,chartX-45);
			var ty=chartY+chartH/2;
			var titleY=canvas.text(tx,ty,C.axisYtext).attr({'font-weight':'bold','font-size':14}).rotate(-90,tx,ty);
			if (C.axisYtitle) C.axisYtitle=titleY; // used by pcaPlot.js to redraw title
		}
		if (C.axisX2text) {
			var titleX2=canvas.text(chartX+chartW/2,chartY-30,C.axisX2text).attr({'font-weight':'bold','font-size':14}); // Math.max(10,chartY-30)
			if (C.axisX2title) C.axisX2title=titleX2;
		}
		if (C.axisY2text) {
			var tx=chartX+chartW+45;
			var ty=chartY+chartH/2;
			var titleY2=canvas.text(tx,ty,C.axisY2text).attr({'font-weight':'bold','font-size':14}).rotate(-90,tx,ty);
			if (C.axisY2title) C.axisY2title=titleY2;
		}


		/***** Drag selection area *****/
		if (C.dragContext) {
			C.dragArea=canvas.rect(0,0,0,0,0)
			.attr({stroke:'#F00',fill:'#600','fill-opacity':0.1})
			.data({startX:0,startY:0,status:'off',chart:C})
			//.click(function(){if (modKeyPressed) {clearDragArea(this)}})
			.dblclick(function() {if (this.data('chart').zoomable) {zoomIn(this.data('chart'))} else {clearDragArea(this)}})
			.hide();

			if (C.zoomable) {
				if (!C.hideZoom) {
					/* Zoom info */
					C.zoomText=canvas.text(chartX,10,'').attr({'font-weight':'bold','font-size':12,'text-anchor':'start'});
				}
				/* Drag area dblclick */
				//C.dragArea.dblclick(function() {zoomIn(this.data('chart'))});
			}
		}

		/***** Thresold lines *****/
		if (C.thresholdLines) {
			for (var i=0; i<C.thresholdLines.length; i++) {
				var thLine=C.thresholdLines[i];
				thLine.pathPos=-50;
				thLine.prevPos=-50;
				if (thLine.axis.match('X')) {
					thLine.path=canvas.path('M'+thLine.pathPos+' '+chartY+' l0 '+chartH)
					.attr({stroke:thLine.color})
					.data({ownerLine:thLine});
					if (thLine.dash) {thLine.path.attr({'stroke-dasharray':thLine.dash});}
					if (thLine.label) {
						thLine.path.hover(
							function(){var th=this.data('ownerLine'); th.popup=drawLabel(canvas,th.pathPos,th.chart.plotArea.attr('y')+th.chart.plotArea.attr('height'),1,th.label,th.color);},
							function(){var th=this.data('ownerLine'); th.popup.remove(); th.popup=null;}
						);
					}
				}
				else { //Y,Y2 axis
					thLine.path=canvas.path('M'+chartX+' '+thLine.pathPos+' l'+chartW+' 0')
					.attr({stroke:thLine.color})
					.data({ownerLine:thLine});
					if (thLine.dash) {thLine.path.attr({'stroke-dasharray':thLine.dash});}
					if (thLine.label) {
						thLine.path.hover(
							function(){var th=this.data('ownerLine'); th.popup=drawLabel(canvas,th.chart.plotArea.attr('x'),th.pathPos,1,th.label,th.color);},
							function(){var th=this.data('ownerLine'); th.popup.remove(); th.popup=null;}
						);
					}
				}
			}
		}
	}

	/***** Data points to SVG points *****/
	if (mainC.dataSets) {
		for (var dsIdx=0; dsIdx < mainC.dataSets.length; dsIdx++) {
			var sizeRule=(mainC.dataSets[dsIdx].params.sizeRule)? mainC.dataSets[dsIdx].params.sizeRule : {min:1,max:20,ratio:2};
			for (var i=0; i < mainC.dataSets[dsIdx].data.length; i++) {
				var dp=mainC.dataSets[dsIdx].data[i];
				//var usedSize=(dp.size)? Math.min(dp.size,20) : 5;
				var usedSize=(dp.size)? Math.max(Math.min(dp.size,sizeRule.max),sizeRule.min) : 5; // max point size
				//var p=canvas.circle(-50,-50,2*Math.sqrt(usedSize/3.14))
				var p=canvas.circle(-50,-50,sizeRule.ratio*Math.sqrt(usedSize/3.14))
				.attr({stroke:'none',fill:mainC.dataSets[dsIdx].params.color,'fill-opacity':mainC.pointOpacity}) //,'fill-opacity':0.65   title:'point '+i,cursor:"crosshair"
				.data({ownerPoint:dp,showLabel:null}) //x:pX,y:pY,
				.hover(function(){setLabelDisplay(this,'on','max')},function(){setLabelDisplay(this,'off')})
				.click(function(){if (modKeyPressed && this.data('ownerPoint').dataSet.params.chart.onPointExclusion) {setPointExclusion(this)} else {setPointSelection(this,'auto','max')}});
				dp.point=p;
				if (dp.valueList) { // multi-point data
					dp.pointList={};
					//var r=Math.round(1.4*Math.sqrt(usedSize/3.14)); // smaller than mater point
					var r=1.4*Math.sqrt(usedSize/3.14); // smaller than mater point
					for (var axis in dp.valueList) {
						dp.pointList[axis]=[];
						for (var j=0; j<dp.valueList[axis].length; j++) {
							dp.pointList[axis].push(canvas.circle(-50,-50,r).attr({stroke:'none',fill:mainC.dataSets[dsIdx].params.color,'fill-opacity':mainC.pointOpacity}).hide()); // hidden
						}
					}
					p.toFront();
				}
			}
		}
	}

    /***** Highlighting legends *****/
	mainC.highlightLegends=canvas.set();

	/***** Chart-specific initialization *****/
	if (mainC.initializeChart) {
		mainC.initializeChart();
	}

	// Looping through displayed charts
	for (var c=0; c<chartList.length; c++) {
		var C=chartList[c];
		/***** Compute default range *****/
		computeDefaultRange(C);

		/***** Plot chart *****/
		plotChart(C);
	}
//,minRangeX,maxRangeX,optRangeX[0],optRangeX[2],minRangeY,maxRangeY,optRangeY[0],optRangeY[2]);
}

function computeDefaultRange(C) { // *public*
	var chartW=C.plotArea.attr('width');
    var chartH=C.plotArea.attr('height');
    var ref=C.chartSettings.reference;
	var minValues={}; var maxValues={};
	if (C.sameScale) { // same scale for both axes at start (PCA)
		var maxRange=0;
		var ranges={};
		for (var axis in C.usedAxes) {
			ranges[axis]=C['maxValue'+axis]-C['minValue'+axis];
			if (maxRange < ranges[axis]) maxRange=ranges[axis];
		}
		for (var axis in ranges) {
			var extraRange=(maxRange-ranges[axis])/2;
			minValues[axis]=C['minValue'+axis]-extraRange;
			maxValues[axis]=C['maxValue'+axis]+extraRange;
		}
	}
	else {
		for (var axis in C.usedAxes) {
			minValues[axis]=C['minValue'+axis];
			maxValues[axis]=C['maxValue'+axis];
		}
	}

	for (var axis in C.usedAxes) {
		var axisSize=(axis.match('X'))? chartW : chartH;
		var optRange=getChartScaleRange(true,minValues[axis],maxValues[axis],axisSize);
		ref[axis]={};
		ref[axis].minRange=(!C['force'+axis+'to0'])? optRange[0] : (C['force'+axis+'to0']==1)? 0 : Math.max(optRange[0],-optRange[3]/5); // 0 or null: auto, 1: 0, 2: ~0
		ref[axis].optMinRange=optRange[1];
		ref[axis].maxRange=optRange[2];
		ref[axis].tickSize=optRange[3];
		ref[axis].pix2valRatio=(ref[axis].maxRange-ref[axis].minRange)/axisSize;
	}

	// set current=reference
    for (var axis in ref) {
		C.chartSettings.current[axis]={};
		for (var prop in ref[axis]) {
			C.chartSettings.current[axis][prop]=ref[axis][prop];
		}
	}
    ref.zoomX=ref.zoomY=C.chartSettings.current.zoomX=C.chartSettings.current.zoomY=1;
}

function plotChart(C) { // *public*
	var mainC=(C.mainChart)? C.mainChart : C;
    var canvas=mainC.canvas;
    var chartX=C.plotArea.attr('x'), chartX0=chartX-0.5;
    var chartY=C.plotArea.attr('y'), chartY0=chartY-1.5;
    var chartW=C.plotArea.attr('width');
    var chartH=C.plotArea.attr('height');
    var settings=C.chartSettings.current;

    /***** Thresold lines *****/
    if (C.thresholdLines) {
		for (var i=0; i<C.thresholdLines.length; i++) {
			moveThresholdLine(C.thresholdLines[i]);
/*
			var thLine=C.thresholdLines[i];
			var axisMin,axisMax;
			if (thLine.axis.match('X')) {
				axisMin=chartX; axisMax=chartX0+chartW;
				thLine.pathPos=chartX0+Math.round((thLine.value-settings[thLine.axis].minRange)/settings[thLine.axis].pix2valRatio);
			}
			else {
				axisMin=chartY; axisMax=chartY0+chartH;
				thLine.pathPos=axisMax-Math.round((thLine.value-settings[thLine.axis].minRange)/settings[thLine.axis].pix2valRatio);
			}
//console.log(thLine.axis+': '+thLine.pathPos+' ('+thLine.prevPos+') ['+axisMin+'-'+axisMax+']');
			if (thLine.pathPos >= axisMin && thLine.pathPos <= axisMax) {
				if (thLine.axis.match('X')) {thLine.path.translate(thLine.pathPos-thLine.prevPos,0);}
				else {thLine.path.translate(0,thLine.pathPos-thLine.prevPos);}
				thLine.path.show();
				thLine.prevPos=thLine.pathPos;
			}
*/
		}
    }

    /***** Axes ticks *****/
    var posX,posY;
	for (var axis in C.usedAxes) {
		// X,X2 ticks
		if (axis.match('X')) {
			if (C['noTicks'+axis]) {
				if (C['noTicks'+axis] != true) {
					posX=chartX0+Math.round(chartW/2);
					posY=(axis=='X')? chartY0+chartH+10 : chartY0-10;
					canvas.text(posX,posY,C.noTicksX).attr({'font-weight':'bold','font-size':12});
				}
			}
			else {
				var sign;
				if (axis=='X') {sign=1; posY=chartY+chartH+1;}
				else {sign=-1; posY=chartY-1;} // posY 2 px below/above chart
				var tickSizeX=settings[axis].tickSize;
				var fixedX=(tickSizeX < 1)? (Math.round(1/tickSizeX)+'').length : (Math.round(tickSizeX)==tickSizeX)? 0 : 1;
				var tickVx=settings[axis].optMinRange;
				while (tickVx <= settings[axis].maxRange) {
					posX=chartX0+Math.round((tickVx-settings[axis].minRange)/settings[axis].pix2valRatio);
					if (posX >= chartX0) {
						C.chartMarks.push(canvas.path('M'+posX+' '+posY+' l0 '+(sign*5)));
						C.chartMarks.push(canvas.text(posX,posY+(sign*10),tickVx.toFixed(fixedX)));
					}
					tickVx+=tickSizeX;
				}
			}
		}
		// Y,Y2 ticks
		if (axis.match('Y')) {
			if (C['noTicks'+axis]) {
				if (C['noTicks'+axis] != true) {
					posX=(axis=='Y')? chartX0-10 : chartX0+chartW+10;
					posY=chartY0+Math.round(chartH/2);
					canvas.text(posX,posY,C.noTicksY).attr({'font-weight':'bold','font-size':12}).rotate(-90,posX,posY);
				}
			}
			else {
				var sign,txtAnch;
				if (axis=='Y') {sign=-1; posX=chartX0; txtAnch='end';}
				else {sign=1; posX=chartX+chartW+1; txtAnch='start';}
				var tickSizeY=settings[axis].tickSize;
				var fixedY=(tickSizeY < 1)? (Math.round(1/tickSizeY)+'').length : (Math.round(tickSizeY)==tickSizeY)? 0 : 1;
				var tickVy=settings[axis].optMinRange;
				while (tickVy <= settings[axis].maxRange) {
					posY=chartY0+chartH-Math.round((tickVy-settings[axis].minRange)/settings[axis].pix2valRatio);
					if (posY <= chartY0+chartH) {
						C.chartMarks.push(canvas.path('M'+posX+' '+posY+' l'+(sign*5)+' 0'));
						C.chartMarks.push(canvas.text(posX+(sign*6),posY,tickVy.toFixed(fixedY)).attr({'text-anchor':txtAnch}));
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
	var existLines=false;
if (mainC.dataSets) {
    var dataSets=mainC.dataSets;
	for (var dsIdx=0; dsIdx < dataSets.length; dsIdx++) {
		//var pathStrg='',startLine=true;
		var dsX=(dataSets[dsIdx].params.axisX)? dataSets[dsIdx].params.axisX : 'X';
		var dsY=(dataSets[dsIdx].params.axisY)? dataSets[dsIdx].params.axisY : 'Y';
		var pathPoints=[];
		//var pointOut=false; var pointIn=false;
		for (var i=0; i < dataSets[dsIdx].data.length; i++) {
			var dp=dataSets[dsIdx].data[i];
			if (dp.subChart && dp.subChart !== C) continue;
			var p=dp.point;
			// point in range
			if (dp.getX() >= settings[dsX].minRange && dp.getX() <= settings[dsX].maxRange && dp.getY() >= settings[dsY].minRange && dp.getY() <= settings[dsY].maxRange) {
				pointIn=true;
				var pX=Math.round((dp.getX()-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0;
				var pY=chartY0+chartH-Math.round((dp.getY()-settings[dsY].minRange)/settings[dsY].pix2valRatio);
				p.attr({cx:pX,cy:pY});
				if (dataSets[dsIdx].params.visible && !dp.isHidden) {
					p.show();
					if ((mainC.connectPoints || dataSets[dsIdx].params.line) && !dp.noLine) {
						//if (pointOut) {
						//	var prevDp=dataSets[dsIdx].data[i-1];
						//	var prevX=Math.round((prevDp.getX()-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0;
						//	var prevY=chartY0+chartH-Math.round((prevDp.getY()-settings[dsY].minRange)/settings[dsY].pix2valRatio);
						//	prevDp.point.attr({cx:prevX,cy:prevY});
						//	pathPoints.push(prevDp.point);
						//	pointOut=false;
						//}
						pathPoints.push(p);
					}
					if (p.data('showLabel')) { // showLabel can be on but not labelObj due to panning
						if (C.dragContext=='pan' && p.data('labelObj')) {
							var tSet=p.data('labelObj');
							tSet.translate(pX-tSet[0].data('x'),pY-tSet[0].data('y'));
							tSet[0].data({x:pX,y:pY});
							tSet.show();
						}
						else {
							emphasizePoint(p,'on',p.data('showLabel'));
						}
					}
					if (dp.pointList) { // multi-point data
						for (var axis in dp.pointList) {
							var pX2,pY2;
							if (axis=='x') {pY2=pY;} else {pX2=pX;}
							for (var j=0; j<dp.pointList[axis].length; j++) {
								if (axis=='x') {pX2=Math.round((dp.valueList.x[j]-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0;}
								else {pY2=chartY0+chartH-Math.round((dp.valueList.y[j]-settings[dsY].minRange)/settings[dsY].pix2valRatio);}
								dp.pointList[axis][j].attr({cx:pX2,cy:pY2});
							}
						}
					}
				}
			}
			// point out of range
			else {
				if (p.data('labelObj') && dataSets[dsIdx].params.visible && C.dragContext=='pan') {
					p.data('labelObj').hide();
				}
				//pointOut=true;
				// Line => compute out of range point & add to path list
				if ((mainC.connectPoints || dataSets[dsIdx].params.line) && !dp.noLine) { //pointIn &&
					var pX=Math.round((dp.getX()-settings[dsX].minRange)/settings[dsX].pix2valRatio)+chartX0;
					var pY=chartY0+chartH-Math.round((dp.getY()-settings[dsY].minRange)/settings[dsY].pix2valRatio);
					p.attr({cx:pX,cy:pY});
					pathPoints.push(dp.point);
					//pointIn=false;
				}
			}
		}
		if (pathPoints.length) {
//console.log(dsIdx+': '+pathPoints.length);
			//var l=canvas.path(pathStrg).attr({stroke:dataSets[dsIdx].params.color,'stroke-width':1});
			dataSets[dsIdx].line=drawPath(C,mainC.canvas,dataSets[dsIdx].params,pathPoints).toBack(); // move behind points
			existLines=true;
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
}

/******************* Datasets visibility ***********************/
function setDataSetDisplay(mainC,dsIdx,visStatus) { // *public*
	mainC.dataSets[dsIdx].params.visible=visStatus;

	if (visStatus) {
		for (var i=0; i<mainC.dataSets[dsIdx].data.length; i++) {
			var dp=mainC.dataSets[dsIdx].data[i];
			var C=(dp.subChart)? dp.subChart : mainC;
			var minX=C.plotArea.attr('x');
			var maxX=minX+C.plotArea.attr('width')-1;
			var minY=C.plotArea.attr('y');
			var maxY=minY+C.plotArea.attr('height')-1;
			var p=dp.point;
			if (p.attr('cx') >= minX && p.attr('cx') <= maxX && p.attr('cy') >= minY && p.attr('cy') <= maxY) {
				p.show();
				//if (p.data('labelObj')) {p.data('labelObj').show()}
				if (p.data('showLabel')) {
					emphasizePoint(p,'on',p.data('showLabel'));
				}
			}
		}
		if (mainC.dataSets[dsIdx].line) {mainC.dataSets[dsIdx].line.show();} //connected points
	}
	else {
		for (var i=0; i<mainC.dataSets[dsIdx].data.length; i++) {
			var dp=mainC.dataSets[dsIdx].data[i];
			dp.point.hide();
			if (dp.point.data('labelObj')) {dp.point.data('labelObj').hide();}
			if (dp.pointList) { // multi-points
				for (var axis in dp.valueList) {
					for (var j=0; j<dp.valueList[axis].length; j++) {
						dp.pointList[axis][j].hide();
					}
				}
			}
		}
		if (mainC.dataSets[dsIdx].line) {mainC.dataSets[dsIdx].line.hide();} //connected points
	}
}

/******************* Point highlighting ***********************/
function addHighlighting(mainC,hName,hColor,pointSet,matchPattern) { // *public* (use ### in matchPattern as point id)
	if (!mainC.allowHighlight) {
		alert('Highlighting is not allowed on chart #'+mainC.chartID);
		return false;
	}
	if (mainC.highlightedPoints[hName]) {
		alert(hName+' is already used.');
		return false;
	}
	mainC.highlightedPoints[hName]={color:hColor,visible:true,dataPoints:[]};
	var dsIdxList={};
	for (var psIdx in pointSet) {
		if (psIdx==-1) { // use all dataSets
			for (var i=0; i < mainC.dataSets.length; i++) {dsIdxList[i]=-1;}
			break;
		}
		else {dsIdxList[psIdx]=psIdx;}
	}
	// Points update
	for (var dsIdx in dsIdxList) {
		var psIdx=dsIdxList[dsIdx];
		P1:for (var i=0; i < pointSet[psIdx].length; i++) {
			var matchRegExp=(matchPattern)? new RegExp(matchPattern.replace('###',pointSet[psIdx][i])) : null;
			for (var j=0; j < mainC.dataSets[dsIdx].data.length; j++) {
				if ((matchPattern && matchRegExp.exec(mainC.dataSets[dsIdx].data[j].externalID)) || (!matchRegExp && mainC.dataSets[dsIdx].data[j].externalID==pointSet[psIdx][i])) {
					mainC.highlightedPoints[hName].dataPoints.push(mainC.dataSets[dsIdx].data[j]);
					mainC.dataSets[dsIdx].data[j].highlightNames.push(hName);
					mainC.dataSets[dsIdx].data[j].isVisble=true; // in case hideUnmatchedPoints is true & 1st active highlight
					var p=mainC.dataSets[dsIdx].data[j].point;
					p.show(); // in case hideUnmatchedPoints is true & 1st active highlight
					if (p.data('showLabel')) {
						if (p.data('labelObj')) {
							p.data('labelObj').remove();
							p.removeData('labelObj');
							displayPointLabel(p,p.data('showLabel')); // label type
						}
					}
					else {p.attr('fill',hColor);}
					if (!matchPattern) continue P1; // more than 1 match per pointSet
				}
			}
		}
	}

	//Update legends on chart itself
	updateHighlightLegends(mainC);

	//DIV update (after points because of point count display)
	updateHighlightList(mainC);

	return true;
}
function deleteHighlighting(mainC,hName) { // *public*
	var hlPoints=mainC.highlightedPoints[hName].dataPoints;
	for (var i=0; i < hlPoints.length; i++) { // dataPoints
		//Removing this hl from list
		var newHlNames=new Array();
		for (var j=0; j < hlPoints[i].highlightNames.length; j++) {
			if (hlPoints[i].highlightNames[j] != hName) {newHlNames.push(hlPoints[i].highlightNames[j]);}
		}
		hlPoints[i].highlightNames=newHlNames;
		var p=hlPoints[i].point;
		if (mainC.hideUnmatchedPoints && newHlNames.length==0) {
			hlPoints[i].isHidden=true;
			p.hide();
		}
		if (p.data('showLabel')) { // point is selected: color=#f00. Do not change it!
			if (p.data('labelObj')) {
				p.data('labelObj').remove();
				p.removeData('labelObj');
				if (!hlPoints[i].isHidden) {displayPointLabel(p,p.data('showLabel'));} // change only label color
			}
		}
		else {
			var newColor=hlPoints[i].dataSet.params.color; // default
			if (newHlNames.length) {
				for (var h=newHlNames.length-1; h>=0; h--) { // pick the newest visible highlight
					if (mainC.highlightedPoints[newHlNames[h]].visible) {
						newColor=mainC.highlightedPoints[newHlNames[h]].color;
						break;
					}
				}
			}
			p.attr('fill',newColor);
		}
	}
	delete mainC.highlightedPoints[hName];

	if (mainC.hideUnmatchedPoints) {
		var noHighlighting=true;
		for (var name in mainC.highlightedPoints) {noHighlighting=false; break;}
		if (noHighlighting) { // reset hideUnmatchedPoints to false
			for (var dsIdx=0; dsIdx < mainC.dataSets.length; dsIdx++) {
				for (var i=0; i < mainC.dataSets[dsIdx].data.length; i++) {
					mainC.dataSets[dsIdx].data[i].isHidden=false;
					if (mainC.dataSets[dsIdx].params.visible) mainC.dataSets[dsIdx].data[i].point.show();
				}
			}
		}
		mainC.hideUnmatchedPoints=false;
	}

	//Update legends on chart itself
	updateHighlightLegends(mainC);

	//DIV update
	updateHighlightList(mainC);

	//Callback function
	if (mainC.updateHighlight) {
		mainC.updateHighlight.callback(hName,'delete');
	}
}
function editHighlighting(mainC,action,hName) { // *public*
	if (action=='edit') {
		document.getElementById(mainC.divID+'_hlightNewName').value=hName;
		document.getElementById(mainC.divID+'_hlightOldName').value=hName;
		document.getElementById(mainC.divID+'_hlightEditDIV').style.display='block';
	}
	else {document.getElementById(mainC.divID+'_hlightEditDIV').style.display='none';}
}
function applyHighlightNameEdition(mainC) { // *public*
	var newName=document.getElementById(mainC.divID+'_hlightNewName').value;
	for (var name in mainC.highlightedPoints) {
		if (newName==name) {
			alert(newName+' is already used.');
			return;
		}
	}
	var oldName=document.getElementById(mainC.divID+'_hlightOldName').value;
	mainC.highlightedPoints[newName]={color:mainC.highlightedPoints[oldName].color,visible:mainC.highlightedPoints[oldName].visible,dataPoints:[]};
	var hlPoints=mainC.highlightedPoints[oldName].dataPoints;
	for (var i=0; i < hlPoints.length; i++) { // dataPoints
		for (var j=0; j < hlPoints[i].highlightNames.length; j++) {
			if (hlPoints[i].highlightNames[j]==oldName) { //Replace oldName with newName
				hlPoints[i].highlightNames[j]=newName;
				break;
			}
		}
		mainC.highlightedPoints[newName].dataPoints.push(hlPoints[i]); // add data point to new
	}
	mainC.highlightedPoints[oldName].dataPoints=[];
	delete mainC.highlightedPoints[oldName];

	//Update legends on chart itself
	updateHighlightLegends(mainC);

	//DIV update (after points because of point count display)
	updateHighlightList(mainC);

	mainC.updateHighlight.callback(oldName,'edit',newName);
}
function updateHighlightLegends(mainC) { // private
	var canvas=mainC.canvas;
	if (!mainC.legendX) {
		mainC.highlightLegends=canvas.set();  // initialize SVG set
		mainC.legendX=canvas.width; // record initial width
	}
	else {mainC.highlightLegends.remove();} // clear previous legend if any
	var posY=25;
	var numAnnot=0;
	for (var hName in mainC.highlightedPoints) {
		if (!mainC.highlightedPoints[hName].visible) {continue;} // hidden highlighting
		numAnnot++;
		if (numAnnot==1) {
			mainC.highlightLegends.push(canvas.text(mainC.legendX,posY,'Legends:').attr({'font-size':12,'font-weight':'bold','text-anchor':'start'}));
			posY+=17;
		}
		mainC.highlightLegends.push(canvas.rect(mainC.legendX+5,posY-5,10,10).attr({fill:mainC.highlightedPoints[hName].color,stroke:mainC.highlightedPoints[hName].color}));
		mainC.highlightLegends.push(canvas.text(mainC.legendX+20,posY,hName).attr({'font-size':12,'text-anchor':'start'}));
		posY+=15;
	}
	/*Readjust chart size if necessary */
	var newWidth=Math.max(mainC.legendX,mainC.highlightLegends.getBBox().x2+15);
	canvas.setSize(newWidth,canvas.height);
	canvas.bottom.attr({width:newWidth}); // adjust background panel too
}
function updateHighlightList(mainC) { //DIV content update
	var hlNameList=new Array();
	for (var name in mainC.highlightedPoints) {hlNameList.push(name);}
	//hlNameList=hlNameList.sort(); // sort is conflicting with point highlightNames order
	var hCode='';
	if (hlNameList.length) {
		hCode='<INPUT type="checkbox" value="1" onclick="setUnmatchedPointsVisibility(cubiojsCharts['+mainC.chartID+'],this.checked)"';
		if (mainC.hideUnmatchedPoints) {hCode+=' checked';}
		hCode+='>Hide unmatched points<BR>';
		for (var i=0; i<hlNameList.length; i++) {
			name=hlNameList[i];
			hCode+='<INPUT type="checkbox" value="'+name+'" onclick="setHighlightVisibility(cubiojsCharts['+mainC.chartID+'],this.value,this.checked)"';
			if (mainC.highlightedPoints[name].visible) hCode+=' checked';
			hCode+='><A href="javascript:selectHighlightedPoints(cubiojsCharts['+mainC.chartID+'],\''+name+'\')" style="font-size:12px;color:'+mainC.highlightedPoints[name].color+'">'+name+' ('+mainC.highlightedPoints[name].dataPoints.length+')</A>&nbsp;';
			if (mainC.updateHighlight && mainC.updateHighlight.editable) {
				hCode+='<INPUT type="button" value="E" style="font-size:10px;font-weight:bold;width:20px" onclick="editHighlighting(cubiojsCharts['+mainC.chartID+'],\'edit\',\''+name+'\')">';
			}
			hCode+='<INPUT type="button" value="X" style="font-size:10px;font-weight:bold;width:20px" onclick="deleteHighlighting(cubiojsCharts['+mainC.chartID+'],\''+name+'\')"><BR>';
		}
		if (mainC.updateHighlight && mainC.updateHighlight.editable) {
			hCode+='<DIV id="'+mainC.divID+'_hlightEditDIV" style="display:none">'; //<FONT style="font-size:12px;font-weight:bold">Name:</FONT>
			hCode+='<INPUT type="text" id="'+mainC.divID+'_hlightNewName" value="" style="width:120px">';
			hCode+='<INPUT type="button" value="V" style="font-weight:bold;font-size:10px;width:20px" onclick="applyHighlightNameEdition(cubiojsCharts['+mainC.chartID+'])"><INPUT type="button" value="X" style="font-size:10px;font-weight:bold;width:20px" onclick="editHighlighting(cubiojsCharts['+mainC.chartID+'],\'cancel\')">';
			hCode+='<INPUT type="hidden" id="'+mainC.divID+'_hlightOldName" value="">';
			hCode+='</DIV>';
		}
	}
	else {hCode='None';}

	document.getElementById(mainC.divID+'_hlight').innerHTML=hCode;
}
function selectHighlightedPoints(mainC,hName) { // *public*
	for (var i=0; i<mainC.highlightedPoints[hName].dataPoints.length; i++) {
		if (!mainC.highlightedPoints[hName].dataPoints[i].isHidden) {
			selectPoint(mainC.highlightedPoints[hName].dataPoints[i].point,'on','min');
		}
	}
}
function setHighlightVisibility(mainC,hName,visStatus) { // *public*
	mainC.highlightedPoints[hName].visible=visStatus;
	var hlPoints=mainC.highlightedPoints[hName].dataPoints;
	for (var i=0; i < hlPoints.length; i++) { // dataPoints
		var dp=hlPoints[i];
		//if (dp.highlightNames[dp.highlightNames.length-1] != hName) {continue;} // no effects on point color
		if (visStatus==false && !dp.isHidden && mainC.hideUnmatchedPoints) {
			var hidePoint=true;
			for (var j=0; j < dp.highlightNames.length; j++) {
				if (dp.highlightNames[j] != hName && mainC.highlightedPoints[dp.highlightNames[j]].visible) {
					hidePoint=false;
					break;
				}
			}
			dp.isHidden=hidePoint;
		}
		var p=dp.point;
		if (p.data('showLabel')) { // point is selected: color=#f00. Do not change it!
			if (dp.isHidden) {
				selectPoint(p,'off');
				//p.hide();
			}
			else if (p.data('labelObj')) {
				p.data('labelObj').remove();
				p.removeData('labelObj');
				displayPointLabel(p,p.data('showLabel')); // change only label color
				//if (!dp.isHidden) displayPointLabel(p,p.data('showLabel')); // change only label color
			}
		}
		else {
			// highlight is the newest for this point (&& other highlights exist) => update point color
			var newColor;
			if (visStatus==true) { // show highlight
				for (var h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
					if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
						newColor=mainC.highlightedPoints[dp.highlightNames[h]].color;
						dp.isHidden=false;
						break;
					}
				}
				//p.show();
			}
			else { // hide highlight
				newColor=dp.dataSet.params.color; // default
				for (var h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
					if (dp.highlightNames[h]==hName) {continue;}
					if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
						newColor=mainC.highlightedPoints[dp.highlightNames[h]].color;
						break;
					}
				}
				//if (dp.isHidden) p.hide();
			}
			p.attr('fill',newColor);
		}
		if (dp.isHidden) {p.hide();} else {p.show();}
	}

	//Update legends on chart itself
	updateHighlightLegends(mainC);
}
function setUnmatchedPointsVisibility(mainC,hideStatus) {
	mainC.hideUnmatchedPoints=hideStatus;
	for (var dsIdx=0; dsIdx < mainC.dataSets.length; dsIdx++) {
		POINT:for (var i=0; i < mainC.dataSets[dsIdx].data.length; i++) {
			var dp=mainC.dataSets[dsIdx].data[i];
			if (hideStatus==true) { // hide unmatched points
				for (var h=0; h < dp.highlightNames.length; h++) {
					if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
						continue POINT;
					}
				}
				dp.isHidden=true;
				if (dp.point.data('showLabel')) selectPoint(dp.point,'off');
				dp.point.hide();
			}
			else {
				dp.isHidden=false;
				if (mainC.dataSets[dsIdx].params.visible) dp.point.show();
			}
		}
	}
}


/******************* Search ***********************/
function searchDataPoints(mainC,searchText,searchResDivID,extSearchChkID,extSearchJobID) { // *public*
	if (extSearchJobID && (!mainC.extSearchJobID || mainC.extSearchJobID != extSearchJobID)) {return;} // ignore this job (another search was launched after it return)
	mainC.extSearchJobID=undefined;
	var searchResDiv=document.getElementById(searchResDivID);
	if (!searchText || searchText.length < 2) {
		searchResDiv.innerHTML='<FONT style="color:#DD0000">Search is string too short!</FONT>';
		return;
	}
	if (extSearchChkID !== null && document.getElementById(extSearchChkID).checked) {
		mainC.extSearchJobID=Date.now();
		searchResDiv.innerHTML='<FONT style="color:#0000FF">Processing...</FONT>';
		mainC.searchable.externalSearch(searchDataPoints,[mainC,searchText,searchResDivID,null,mainC.extSearchJobID],1); // (local search function,array of function arguments,index of search text in array). To be recalled by external search function
		return;
	}
//	var minVx,maxVx,minVy,maxVy;
	var matchList=[];
	var matchExp=new RegExp(searchText,"i");
	var neverMatched=true, okProceed100=false, okProceed1000=false;
	var newZoom={};
	for (var dsIdx=0; dsIdx<mainC.dataSets.length; dsIdx++) {
		if (!mainC.dataSets[dsIdx].params.visible) continue;
		for (var i=0; i < mainC.dataSets[dsIdx].data.length; i++) {
			var dp=mainC.dataSets[dsIdx].data[i];
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
				if (dp.point.attr('cx') < 0) { // point if off chart ZoomOut required
					var chartName=(dp.subChart)? dp.subChart.name : 'main';
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
	for (var chName in newZoom) {
		if (chName=='main') {zoomOut(mainC,true);}
		else {zoomOut(mainC.subChart[chName],true);}
	}
	for (var i=0; i < matchList.length; i++) {
		//emphasizePoint(matchList[i].point,'on','min');
		selectPoint(matchList[i].point,'on','min');
	}
	searchResDiv.innerHTML=matchList.length+' match(es) found!';
}

/******************* Selection management ***********************/
function manageSelection(mainC,action) { // *public*
	var chartList=[];
	if (mainC.subChart) {
		for (var c=0; c<mainC.activeCharts.length; c++) {
			chartList.push(mainC.subChart[mainC.activeCharts[c]]);
		}
	}
	else {chartList.push(mainC);}

	// Looping through displayed charts
	for (var c=0; c<chartList.length; c++) {
		var C=chartList[c];
		if (!C.dragArea || C.dragArea.data('status')=='off') continue; //return;
		var minPx=C.dragArea.attr('x');
		var maxPx=minPx+C.dragArea.attr('width')-1;
		var minPy=C.dragArea.attr('y');
		var maxPy=minPy+C.dragArea.attr('height')-1;
		for (var dsIdx=0; dsIdx<mainC.dataSets.length; dsIdx++) {
			if (!mainC.dataSets[dsIdx].params.visible) continue;
			for (var i=0; i < mainC.dataSets[dsIdx].data.length; i++) {
				var dp=mainC.dataSets[dsIdx].data[i];
				if (dp.isHidden || (dp.subChart && dp.subChart !== C)) continue;
				var p=dp.point;
				if (p.attr('cx')>=minPx && p.attr('cx')<=maxPx && p.attr('cy')>=minPy && p.attr('cy')<=maxPy) {
					selectPoint(p,action,'min');
					//emphasizePoint(p,action,'min');
				}
			}
		}
	}
}
function listSelectedPoints(C,actionIndex) { // *public*
	var selectedPoints=new Object();
	var count=0;
	for (var dsIdx=0; dsIdx<C.dataSets.length; dsIdx++) {
		if (!C.dataSets[dsIdx].params.visible) continue;
		selectedPoints[dsIdx]=new Array();
		for (var i=0; i < C.dataSets[dsIdx].data.length; i++) {
			var dp=C.dataSets[dsIdx].data[i];
			if (dp.point.attr('cx') > 0 && dp.point.data('labelObj')) { // visible & selected
				selectedPoints[dsIdx].push(dp.externalID);
				count++;
			}
		}
		if (selectedPoints[dsIdx].length==0) {delete selectedPoints[dsIdx];}
	}
	if (count==0) {alert('No selected points found in displayed range.');}
	else {
		var thresholds=new Array();
		for (var i=0; i<C.thresholdLines.length;i++) {
			thresholds.push(C.thresholdLines[i].initValue);
		}
		if (actionIndex < 0) {
			C.pointOnList(selectedPoints,thresholds);
		}
		else {
			C.pointOnList[actionIndex][1](selectedPoints,thresholds);
		}
	}
}

/******************* Thresholds management ***********************/
function selectThreshold(C,thIdx,thresBoxID) { // *public*
	document.getElementById(thresBoxID).value=(thIdx>=0)? C.thresholdLines[thIdx].initValue : '';
}
function updateThreshold(mainC,thresSelID,thresBoxID) { // *public*
	var thIdx=document.getElementById(thresSelID).value;
	var newValue=document.getElementById(thresBoxID).value;
	if (!isNumber(newValue)) {alert(newValue+' is not a valid number!'); return;}
	if (thIdx<0 || !newValue) return;
	mainC.thresholdLines[thIdx].setValues(newValue);
	mainC.thresholdLines[thIdx].path.hide();
	//plotChart(mainC.thresholdLines[thIdx].chart);
	moveThresholdLine(mainC.thresholdLines[thIdx]);

}
function moveThresholdLine(thLine) {
	var mainC=(thLine.chart.mainChart)? thLine.chart.mainChart : thLine.chart;
    var canvas=mainC.canvas;
    var chartX=thLine.chart.plotArea.attr('x'), chartX0=chartX-0.5;
    var chartY=thLine.chart.plotArea.attr('y'), chartY0=chartY-1.5;
    var chartW=thLine.chart.plotArea.attr('width');
    var chartH=thLine.chart.plotArea.attr('height');
    var settings=thLine.chart.chartSettings.current;
	var axisMin,axisMax;
	if (thLine.axis.match('X')) {
		axisMin=chartX; axisMax=chartX0+chartW;
		thLine.pathPos=chartX0+Math.round((thLine.value-settings[thLine.axis].minRange)/settings[thLine.axis].pix2valRatio);
	}
	else {
		axisMin=chartY; axisMax=chartY0+chartH;
		thLine.pathPos=axisMax-Math.round((thLine.value-settings[thLine.axis].minRange)/settings[thLine.axis].pix2valRatio);
	}
	if (thLine.pathPos >= axisMin && thLine.pathPos <= axisMax) {
		if (thLine.axis.match('X')) {thLine.path.translate(thLine.pathPos-thLine.prevPos,0);}
		else {thLine.path.translate(0,thLine.pathPos-thLine.prevPos);}
		thLine.path.show();
		thLine.prevPos=thLine.pathPos;
	}
}
function isNumber(n) { // *public*
  return !isNaN(parseFloat(n)) && isFinite(n);
}

/******************* Chart zooming ***********************/
function zoomIn(C) {
var mainC=(C.mainChart)? C.mainChart : C;
    var chartX0=C.plotArea.attr('x')-1;
    var chartY0=C.plotArea.attr('y')-1;
    var chartW=C.plotArea.attr('width');
    var chartH=C.plotArea.attr('height');
	var minPx,maxPx,minPy,maxPy;
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
    var ref=C.chartSettings.reference;
    var cur=C.chartSettings.current;

	for (var axis in C.usedAxes) {
		var minV,maxV,axisSize;
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
		var newRange=getChartScaleRange(false,minV,maxV,axisSize);

		cur[axis].minRange=minV;
		cur[axis].optMinRange=newRange[1];
		cur[axis].maxRange=maxV;
		cur[axis].tickSize=newRange[3];
		cur[axis].pix2valRatio=(cur[axis].maxRange-cur[axis].minRange)/axisSize;

		var axisZoom=(axis.match('X'))? 'zoomX' : 'zoomY';
		cur[axisZoom]=ref[axis].pix2valRatio/cur[axis].pix2valRatio;
	}

    /*** Clear & plot chart ***/
    clearChart(C);
    plotChart(C);
}

function zoomOut(C,zoomTo1) { // *public*
    var cur=C.chartSettings.current;
    if (cur.zoomX==1 && cur.zoomY==1) {return;}
    var ref=C.chartSettings.reference;

	for (var axis in C.usedAxes) {
		var minV,maxV,axisSize,axisZoom;
		if (modKeyPressed || zoomTo1) { // set zooms to 1
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
}


/******************* Chart panning ***********************/
function panChart(C) {
	var chartW=C.plotArea.attr('width');
	var chartH=C.plotArea.attr('height');
	var ref=C.chartSettings.reference;
    var cur=C.chartSettings.current;

	for (var axis in C.usedAxes) {
		var dV,axisSize;
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
		var newRange=getChartScaleRange(false,cur[axis].minRange,cur[axis].maxRange,axisSize);
		cur[axis].optMinRange=newRange[1];
		cur[axis].tickSize=newRange[3];
	}

    /*** Clear & plot chart ***/
    clearChart(C);
    plotChart(C);
}

function clearChart(C) { // *public*
    clearDragArea(C.dragArea);

if (C.thresholdLines) {
    for (var i=0; i<C.thresholdLines.length; i++) {
		if (C.thresholdLines[i].path) {
			C.thresholdLines[i].path.hide();
			/*
			C.thresholdLines[i].path.remove();
			C.thresholdLines[i].path=null;
			C.thresholdLines[i].pathPos=null;
			if (C.thresholdLines[i].popup) {
				C.thresholdLines[i].popup.remove();
				C.thresholdLines[i].popup=null;
			}
			*/
		}
    }
}
	for (var i=0; i<C.chartMarks.length; i++) {
		C.chartMarks[i].remove();
	}
    C.chartMarks.length=0;

//var maxDx=C.plotArea.attr('width')+50;
//var maxDy=C.plotArea.attr('height')+50;

	var mainC=(C.mainChart)? C.mainChart : C;
if (mainC.dataSets) {
    for (var dsIdx=0; dsIdx<mainC.dataSets.length; dsIdx++) {
		for (var i=0; i<mainC.dataSets[dsIdx].data.length; i++) {
			var dp=mainC.dataSets[dsIdx].data[i];
			if (dp.subChart && dp.subChart !== C) continue;
			var p=dp.point;
			p.hide();
			p.attr({cx:-100,cy:-100});
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
		if (mainC.dataSets[dsIdx].line) mainC.dataSets[dsIdx].line.remove();
    }
}
	if (C.clearChart) {
		C.clearChart();
	}
}

function getChartScaleRange(optimize,minV,maxV,dimSize) { // dimSize is optional  // *public*
	var deltaV=maxV-minV;
	if (deltaV==0) { // minV=maxV => 1 datapoint
		minV-=0.05*minV;
		maxV+=0.05*maxV;
		deltaV=maxV-minV;
	}
	else if (optimize) {
		minV-=(0.05*deltaV); // - 5% range
		maxV=(maxV*1)+(0.05*deltaV); // + 5% range (force to number)
		deltaV*=1.1;
	}
//console.log('Min='+minV+', Max='+maxV+', Delta='+deltaV);
	//Move delta within ]1-10] range
	var shiftFactor=1;
	if (deltaV <= 1) {
		//while (deltaV*shiftFactor < 1) {shiftFactor*=10;}
		var v=Math.floor(1/deltaV)+'';
//console.log(1/deltaV+' => '+v);
		shiftFactor=Math.pow(10,v.length);
	}
	else if (deltaV > 10) {
		//while (deltaV*shiftFactor > 10) {shiftFactor/=10;}
		var v=Math.floor(deltaV - 1)+''; // 100 => 99 (length=2 !3)
//console.log(deltaV+' => '+v);
		shiftFactor=Math.pow(10,-(v.length - 1));
	}
	var optMinV=Math.floor(minV*shiftFactor)/shiftFactor;
	//var maxV=maxV+(0.05*deltaV); // + 5% range
	var scaledRange=(maxV-optMinV)*shiftFactor; //deltaV*shiftFactor;
	var scaledTick;
	if (dimSize && dimSize < 250) { // 5 ticks
		scaledTick=(scaledRange > 5)? 2 : (scaledRange > 2.5)? 1 : (scaledRange > 1.1)? 0.5 : 0.2;
	}
	else { // 10 ticks
		scaledTick=(scaledRange > 5)? 1 : (scaledRange > 2.5)? 0.5 : (scaledRange > 1.1)? 0.2 : 0.1;
	}
	var tickSize=scaledTick/shiftFactor;
//console.log('D='+deltaV+', F='+shiftFactor+', minV='+minV+', optMinV='+optMinV+', tick='+tickSize);
	return [minV,optMinV,maxV,tickSize];
}


/******************* Point selection & highlight ***********************/
function setPointSelection(p,action,type) { // point click event // *public*
	var mainC=p.data('ownerPoint').dataSet.params.chart;
	if (mainC.selectAllDataSets) {
		var linkedPoints=mainC.dataSetLinks[p.data('ownerPoint').externalID];
		var refAction=(action=='auto' && p.data('showLabel'))? 'off' : (action=='auto')? 'on' : action; // set action based on reference point
		for (var i=0; i < linkedPoints.length; i++) {
			selectPoint(linkedPoints[i].point,refAction,type);
		}
	}
	else {selectPoint(p,action,type);}
}
function selectPoint(p,action,type) {
    if (action=='auto') { // invert selection
		if (p.data('showLabel')) { // already showing label => set OFF
			p.data('showLabel',null);
			emphasizePoint(p,'off');
		}
		else { // set ON
			p.data('showLabel',type);
			emphasizePoint(p,'on',type);
		}
	}
	else if (action=='on') { // set ON
		if (!p.data('showLabel')) {
			p.data('showLabel',type);
			emphasizePoint(p,'on',type);
		}
	}
	else { // OFF
		p.data('showLabel',null);
		emphasizePoint(p,'off');
	}
    //alert('point ['+p.id+'] '+p.data('index')+' ('+p.data('x')+','+p.data('y')+')');
}
function setPointExclusion(p) { // *public*
	selectPoint(p,'off');
	var dp=p.data('ownerPoint');
	dp.noLine=(dp.noLine)? false : true;
	var dSet=dp.dataSet;
	dSet.line.remove();
	//var pathStrg='',startLine=true;
	var pathPoints=[];
	var dpIdx;
	for (var i=0; i < dSet.data.length; i++) {
		var dpi=dSet.data[i];
		var pi=dpi.point;
		if (!dpi.noLine) {
			pathPoints.push(pi);
			//if (startLine==true) {
			//	pathStrg='M'+pi.attr('cx')+','+pi.attr('cy')+' R'; // curve to
			//	startLine=false;
			//}
			//else {pathStrg+=' '+pi.attr('cx')+','+pi.attr('cy');}
		}
		if (dpi===dp) {dpIdx=i;}
	}
	var mainC=dSet.params.chart;
	if (mainC.connectPoints || dSet.params.line) {
		//dSet.line=C.canvas.path(pathStrg).attr({stroke:dSet.params.color,'stroke-width':1}).toBack();
		var C=(dp.chart)? dp.chart : mainC;
		dSet.line=drawPath(C,mainC.canvas,dSet.params,pathPoints).toBack();
		C.plotArea.toBack();
		C.backPanel.toBack();
	}
//console.log(dSet.params.index+','+dpIdx);
	if (mainC.onPointExclusion) {mainC.onPointExclusion(dSet.params.index,dpIdx,dp.noLine);}
}
function setLabelDisplay(p,action,type) { // point hover event  // *public*
	var C=p.data('ownerPoint').dataSet.params.chart;
	if (C.dataSets.length > 1 && C.dataSetLinks && p.data('ownerPoint').externalID && C.dataSetLinks[p.data('ownerPoint').externalID]) { // not defined for all charts
		var linkedPoints=C.dataSetLinks[p.data('ownerPoint').externalID];
		for (var i=0; i < linkedPoints.length; i++) {
			var usedType=(C.selectAllDataSets || linkedPoints[i]===p.data('ownerPoint'))? type : 'none';
			emphasizePoint(linkedPoints[i].point,action,usedType);
		}
	}
	else {emphasizePoint(p,action,type);}
}
function emphasizePoint(p,action,type) { // point hover & click events
    var pColor;
	var dp=p.data('ownerPoint');
	var mainC=dp.dataSet.params.chart;
	var C=(dp.subChart)? dp.subChart : mainC;
	if (action=='on') {
		pColor=(mainC.highlightColor)? mainC.highlightColor : '#F00';
		//if (p.data('labelObj')) return;
		/* Checking if point is visible (hide if not) */
		if (type != 'none') {
			var minX=C.plotArea.attr('x');
			var maxX=minX+C.plotArea.attr('width')-1;
			var minY=C.plotArea.attr('y');
			var maxY=minY+C.plotArea.attr('height')-1;
			if (dp.dataSet.params.visible && p.attr('cx') >= minX && p.attr('cx') <= maxX && p.attr('cy') >= minY && p.attr('cy') <= maxY) {
				if (p.data('labelObj')) {p.data('labelObj').show();}
				else {displayPointLabel(p,type);}
			}
			if (dp.pointList) {
				for (var axis in dp.valueList) {
					for (var j=0; j<dp.valueList[axis].length; j++) {
						dp.pointList[axis][j].show();
					}
				}
			}
		}
	}

	else { // off
		if (p.data('showLabel')) { // needed when hovering out of point
			pColor=(mainC.highlightColor)? mainC.highlightColor : '#F00';
		}
		else {
			//pColor=(dp.highlightNames.length)? mainC.highlightedPoints[dp.highlightNames[dp.highlightNames.length-1]].color : dp.dataSet.params.color;
			pColor=dp.dataSet.params.color; // default
			if (dp.highlightNames) {
				for (var h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
					if (mainC.highlightedPoints[dp.highlightNames[h]].visible) {
						pColor=mainC.highlightedPoints[dp.highlightNames[h]].color;
						break;
					}
				}
			}
if (type != 'none') {
			if (p.data('labelObj')) { // !!!! Not sure this can happen???? (13/02/13)
				p.data('labelObj').remove();
				p.removeData('labelObj');
			}
			if (dp.pointList) {
				for (var axis in dp.valueList) {
					for (var j=0; j<dp.valueList[axis].length; j++) {
						dp.pointList[axis][j].hide();
					}
				}
			}
}
		}
	}
	p.attr('fill',pColor);
}
function displayPointLabel(p,type) { // !!! also in peptidePlot.js !!!
	var dp=p.data('ownerPoint');
	var C=dp.dataSet.params.chart;
	var text=(C.customPointLabel)? C.customPointLabel(dp,type) : dp.info(type); //user-defined text
	//var lColor=(dp.highlightNames.length)? C.highlightedPoints[dp.highlightNames[dp.highlightNames.length-1]].color : dp.dataSet.params.color;
	var lColor=dp.dataSet.params.color;// default
	if (dp.highlightNames && dp.highlightNames.length) {
		for (var h=dp.highlightNames.length-1; h>=0; h--) { // pick the newest visible highlight
			if (C.highlightedPoints[dp.highlightNames[h]].visible) {
				lColor=C.highlightedPoints[dp.highlightNames[h]].color;
				break;
			}
		}
	}

	var tSet=drawLabel(p.paper,p.attr('cx'),p.attr('cy'),p.attr('r')+1,text,lColor);
	if (dp.dataSet.params.chart.pointOnclick) {
		tSet.click(function(){dp.dataSet.params.chart.pointOnclick(dp.dataSet.params.index,dp.label,dp.externalID)})
		.attr({cursor:'pointer'});
	}
	tSet[0].data({x:p.attr('cx'),y:p.attr('cy')});
//tSet.translate(0,50);
//console.log('SET: x='+tSet[0].data('x')+', y='+tSet[0].data('y'));
	p.data('labelObj',tSet);
}
function drawLabel(canvas,x,y,d,text,lColor) { // !!! also in peptidePlot.js !!! // *public*
//console.log(lColor);
	var shift=15;
	var t=canvas.text(x+d+shift,y-d-shift,text).attr({'font-size':10,'text-anchor':'start',fill:lColor}); //,'font-weight':'bold','fill-opacity':0.6
	var tBox=t.getBBox();
	var tx=tBox.x;
	var ty=tBox.y;
	var tw=tBox.width;
	var th=tBox.height;
	//var tb; //Bubble around text
	//var dy=Math.round((th-6)/2);
	var th05=Math.round(th/2);
	var sx;
	if (tx+tw > canvas.width) {
		//t.attr({x:canvas.width-tBox.width});
		tx=tx-2*(d+shift)-tw; //tBox.width;
		t.attr({x:tx});
		//tb=canvas.path('M'+(tx-2)+' '+(ty-1)+' l'+(tw+4)+' 0 l0 '+dy+' l10 3 l-10 3 l0 '+dy+' l'+(-tw-4)+' 0 Z');
		sx=tx+tw+3;
	}
	else {
		//tb=canvas.path('M'+(tx-2)+' '+(ty-1)+' l'+(tw+4)+' 0 l0 '+(th+2)+' l'+(-tw-4)+' 0 l0 '+(-dy)+' l-10 -3 l10 -3 Z');
		sx=tx-1;
	}
	if (ty < 0) {
	    t.attr({y:2+th05});
	    ty=t.getBBox().y;
	}
	var sy=ty+th05;
	var ex=x+Math.round(d*(sx-x)/(d+shift));
	var ey=y+Math.round(d*(sy-y)/(d+shift));
	var tl=canvas.path('M'+sx+' '+sy+' L'+ex+' '+ey).attr({stroke:lColor});
	var tb=canvas.rect(Math.round(tx)-1.5,Math.round(ty)-1.5,tw+4,th+2,4).attr({stroke:lColor,fill:'#fff','fill-opacity':0.7});
	t.toFront();

	return canvas.set(tb,tl,t);
}

/********** Path drawing **********/
function drawPath(C,canvas,params,pathPoints) { // public
	var type='curve';
	var color=params.color;
	var pattern;
	var width=1;
	var opacity=1;
	if (params.line) {
		if (params.line.type) type=params.line.type;
		if (params.line.color) color=params.line.color;
		if (params.line.pattern) pattern=params.line.pattern;
		if (params.line.width) width=params.line.width;
		if (params.line.opacity) opacity=params.line.opacity;
	}

	var pathStrg='';
	if (type=='function') {

	}
	else {
		if (type != 'line' && pathPoints.length <= 2) type='line';
		for (var i=0; i<pathPoints.length; i++) {
			if (i==0) {
				pathStrg='M'+pathPoints[i].attr('cx')+' '+pathPoints[i].attr('cy');
				pathStrg+=(type=='curve')? ' R' : ' L';
			}
			else {pathStrg+=' '+pathPoints[i].attr('cx')+' '+pathPoints[i].attr('cy');}
		}
	}
	var clipStrg=C.plotArea.attr('x')+','+C.plotArea.attr('y')+','+C.plotArea.attr('width')+','+C.plotArea.attr('height');
	var l=canvas.path(pathStrg).attr({stroke:color,'stroke-width':width,'stroke-opacity':opacity,'clip-rect':clipStrg});
	if (pattern) {l.attr({'stroke-dasharray':pattern});}

for (var i=0; i<pathPoints.length; i++) {pathPoints[i].toFront();}
	return l;
}

/******************* Mouse drag management *********************/
function setDragging(C,x,y,e) { // x,y are absolute to document // *public*
    //var mouseInChart=getMousePositionInChart(e);
	var mousePos=getMousePositionInElement(e);
//debug('ctrl='+modKeyPressed+' : '+x+'('+mousePos[0]+'), '+y+'('+mousePos[1]+')');
	var mouseX=mousePos[0]+C.plotArea.attr('x')-1,mouseY=mousePos[1]+C.plotArea.attr('y')-1;
	/*
	if (e.offsetX) { // IE
		if (e.layerX) { // IE 9+
			mouseX=e.offsetX;
			mouseY=e.offsetY;
		}
		else { // IE 5-8
			mouseX=e.offsetX+C.plotArea.attr('x')-1;
			mouseY=e.offsetY+C.plotArea.attr('y')-1;
		}
	}
	else {
		mouseX=e.layerX;
		mouseY=e.layerY;
	}
	*/
    if (modKeyPressed) {
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
		//setDragArea(C,mouseInChart[0],mouseInChart[1]);
		setDragArea(C,mouseX,mouseY);
    }
}
function extendDragging(C,dx,dy,x,y,e) { // *public*
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
    else {extendDragArea(C,dx,dy);}
}
function endDragging(C) { // *public*
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
		endDragArea(C);
//		C.dragContext=null;
    }
    C.plotArea.attr({cursor:'auto'});
}
/******** Drag area management *******/
function setDragArea(C,mX,mY) { // position in chart
var mainC=(C.mainChart)? C.mainChart : C;
var x,y,w,h;
if (C.noZoomX && !mainC.dataSets) {x=C.plotArea.attr('x'); w=C.plotArea.attr('width');} //
else {x=mX; w=0;}
if (C.noZoomY && !mainC.dataSets) {y=C.plotArea.attr('y'); h=C.plotArea.attr('height');} //
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
function extendDragArea(C,dx,dy) {
var mainC=(C.mainChart)? C.mainChart : C;
    var dragArea=C.dragArea;
    var plotArea=C.plotArea;
    var chartX0=plotArea.attr('x')-1;
    var chartY0=plotArea.attr('y')-1;
	if (!C.noZoomX || mainC.dataSets) {
    if (dx > 0) {
		dragArea.attr({width:dx});
		if (dragArea.attr('x')+dragArea.attr('width') > plotArea.attr('width')+chartX0) {dragArea.attr({width:plotArea.attr('width')-dragArea.attr('x')+chartX0})};
    }
    else {
		dragArea.attr({x:dragArea.data('startX')+dx,width:-dx});
		if (dragArea.attr('x') < plotArea.attr('x')) {
			var extra=plotArea.attr('x') - dragArea.attr('x');
			dragArea.attr({x:plotArea.attr('x'),width:dragArea.attr('width')-extra});
		}
    }
}
if (!C.noZoomY || mainC.dataSets) {
    if (dy > 0) {
	    dragArea.attr({height:dy});
	    if (dragArea.attr('y')+dragArea.attr('height') > plotArea.attr('height')+chartY0) {dragArea.attr({height:plotArea.attr('height')-dragArea.attr('y')+chartY0})};
    }
    else {
		dragArea.attr({y:dragArea.data('startY')+dy,height:-dy});
		if (dragArea.attr('y') < plotArea.attr('y')) {
			var extra=plotArea.attr('y') - dragArea.attr('y');
			dragArea.attr({y:plotArea.attr('y'),height:dragArea.attr('height')-extra});
		}
    }
}
//console.log('EXT:'+dragArea.attr('x')+','+dragArea.attr('y')+','+dragArea.attr('width')+','+dragArea.attr('height'));
}
function endDragArea(C) {
var mainC=(C.mainChart)? C.mainChart : C;
if ((!mainC.dataSets && ((C.noZoomX && C.dragArea.attr('height') > 5) || (C.noZoomY && C.dragArea.attr('width') > 5))) || ((mainC.dataSets || !C.noZoomX && !C.noZoomY) && C.dragArea.attr('width')*C.dragArea.attr('height') > 25)) { // big enough => end extension
		C.dragArea.data({status:'on'});
if (C.autoZoom) {
	zoomIn(C);
}
    }
    else { // too small => hide
		clearDragArea(C.dragArea);
    }
}

function clearDragArea(dragArea) {
	dragArea.hide().data({status:'off'}).attr({width:0,height:0});
}

/***************** Mouse and Keyboard events ********************/
document.onkeydown=function(e){
	var kbEv=window.event || e;
	if (kbEv.keyCode==16 || kbEv.keyCode==17 || kbEv.keyCode==18) { // ctrl or alt
		modKeyPressed=true;
	}
}
document.onkeyup=function(e){
	var kbEv=window.event || e;
	if (kbEv.keyCode==16 || kbEv.keyCode==17 || kbEv.keyCode==18) { // ctrl or alt
		modKeyPressed=false;
	}
}

//function getMousePositionInChart(e) {
//	var mousePos=e.layerX? [e.layerX,e.layerY] : [e.offsetX+chartX0,e.offsetY+chartY0];
//	return mousePos;
//}
function getMousePositionInElement(e) {
	var target = e.target || e.srcElement,
	style = target.currentStyle || window.getComputedStyle(target, null),
	borderLeftWidth = parseInt(style['borderLeftWidth'], 10),
	borderTopWidth = parseInt(style['borderTopWidth'], 10),
	rect = target.getBoundingClientRect();
	return [e.clientX - borderLeftWidth - rect.left,e.clientY - borderTopWidth - rect.top];
}

/************* String processing *************/
function shortenText(myText,maxLength) {
	var textLength=myText.length;
	if (textLength <= maxLength) return myText;
	var subTextSize=Math.max(2,Math.floor((maxLength-3)/2));
	return myText.substring(0,subTextSize)+'...'+myText.substring(textLength-subTextSize);
}

/************* Chart Export *************/
function exportSVGtoImg (svgDivId,imgName,exportScript,format) {
	if (!format) format='png';

	//put html svg to svgHTML
	//svgFix make svg standard
	var svgHTML = document.getElementById(svgDivId).innerHTML;
	var data;
	if (format=='png') { // specific processing for png using canvas element
		//create canvas element
		var canvas = document.createElement('canvas');
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
					exportDataStringImage();
				}
			}
		);
	}
	else if (format=='svg') {
		data = svgfix(svgHTML);
		exportDataStringImage();
	}
	function exportDataStringImage() {
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
}


/*
####>Revision history<####
# 1.5.0 [FEATURE] Search: Optional external user-provided pre-search (PP 04/10/19)
# 1.4.5 Double click removes drag area on non-zoomable charts (PP 15/05/19)
# 1.4.4 Displays both PNG and SVG image export options if format is not specified by user (PP 06/06/18)
# 1.4.3 Improved chart.showDataSets detection & chart-specific form elements (PP 20/11/16)
# 1.4.2 Bug fix in new pointOnList behavior (PP 10/11/16)
# 1.4.1 pointOnList now also accepts an array of [action,function] array (PP 03/11/16)
# 1.4.0 Drawing legends for highlights (PP 23/10/16)
# 1.3.8 Better handling of single data poin in getChartScaleRange (PP 27/05/16)
# 1.3.7 Added chart registering (PP 20/02/16)
# 1.3.6 max-height on Datasets FIELDSET (PP 28/10/15)
# 1.3.5 Handle dataSet.params.sizeRule for scatter plots (PP 21/10/15)
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
# 1.1.0 Better check on point linkage across dataSets during mouseover (PP 29/01/13)
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
