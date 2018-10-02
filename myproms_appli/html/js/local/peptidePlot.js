/*
################################################################################
# peptidePlot.js             1.1.4                                             #
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
                  Peptide Plot object
******************************************************************/
function peptidePlot(plotData) {
	var PP=this;
    var TOP_SPACE=20.5;
    var RIGHT_SPACE=20.5;
    var XAXIS_SPACE=50.5;
    var YAXIS_SPACE=70.5;
	//this.pointOnclick=plotData.pointOnclick;
	//this.pointOnList=plotData.pointOnList;
    var div=document.getElementById(plotData.div);
    var chartW=(plotData.width)? plotData.width : 400;
    var chartH=(plotData.height)? plotData.height : 400;
	var chartX0=YAXIS_SPACE;
    var chartX=YAXIS_SPACE+1;
	var chartY0=TOP_SPACE;
    var chartY=TOP_SPACE+1;
    var canvasW=YAXIS_SPACE+chartW+RIGHT_SPACE; // can be modified below
    var canvasH=TOP_SPACE+chartH+XAXIS_SPACE;
	var chartYLabel=plotData.valueAxisLabel || 'Values';
	var valueLabel=plotData.valueLabel || 'Value';
	var convertValue=plotData.convertValue || function(v) {return v}; // PP. because needed in peptide object
	var convertValueDisplayed=plotData.convertValueDisplayed || function(v) {return v.toPrecision(3)};
	var peptideProperties=plotData.peptideProperties || [];
	var pepIdPropertyIndex;
	for (var i=0; i < peptideProperties.length; i++) {
		if (peptideProperties[i]=='id') {
			pepIdPropertyIndex=i;
			break;
		}
	}

    var minValueY=null,maxValueY=null;

	/******* Protein data ******/
	var protein={
		name: plotData.protein.name || 'Protein',
		length: plotData.protein.length || 0,
		value: plotData.protein.value || null, // can be undefined
		valueY: (plotData.protein.value)? convertValue(plotData.protein.value) : null,
		relativeThresholds:	plotData.protein.relativeThresholds || null,
		thresholdsY: null,
		popup: null
    };
	if (protein.value != null) {
		minValueY=protein.valueY;
		maxValueY=protein.valueY;
		if (protein.relativeThresholds) {
			protein.thresholdsY=[];
			for (var i=0; i<protein.relativeThresholds.length; i++) {
				var cValue=convertValue(protein.relativeThresholds[i] * protein.value);
				protein.thresholdsY.push(cValue);
				minValueY=Math.min(minValueY,cValue);
				maxValueY=Math.max(maxValueY,cValue);
			}
		}
	}

	/******* Reference line (eg. to show fold change=1) ******/
	var reference={};
	if (plotData.reference && plotData.reference.value != undefined) {
		reference.label=plotData.reference.label || 'Reference';
		reference.value=plotData.reference.value;
		reference.valueY=convertValue(plotData.reference.value);
		minValueY=Math.min(minValueY,reference.valueY);
		maxValueY=Math.max(maxValueY,reference.valueY);
	}

	if (protein.value || reference.value) canvasW+=50; // increase right-side space to write protein value

    /********************* Peptide import *********************/
    var peptideSet=[]; // list of peptides
    var lastAAcovered=1;
	var peptideColor=plotData.peptideColor || {};
	if (peptideColor.map && peptideColor.map.toLowerCase()=='value') { // matches peptide primary value
		peptideColor.map=valueLabel;
	}
	var noColorRange=true;
	if (peptideColor.range) { // incompatible with peptideColor.type=discrete and peptideColor.list
		noColorRange=false;
		if (peptideColor.map==valueLabel) { // matches peptide primary value => convert range
			peptideColor.range[0]=convertValue(peptideColor.range[0]);
			peptideColor.range[1]=convertValue(peptideColor.range[1]);
			peptideColor.type='continuous'; // just in case
			peptideColor.list=null;  // just in case
		}
	}
	var colorList=['#0000FF','#4AA02C','#000000','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C','#E18B6B'];
	var mappedPropertyIdx=null;
	var mappedValues={}; // only used is color mapping on discrete values
	var mappingError=false;

	function peptide(data) { /* Peptide object */
		this.startPos=(typeof(data[0])=='number')? [data[0]] : data[0]; // single position OR array of positions
		this.sequence=data[1];
		this.length=data[1].length;
		this.varMod=data[2]; // can be null
		this.charge=data[3]; // can be null
		this.value=data[4];
		this.valueY=convertValue(data[4]);
		this.excluded=(data[5])? true : false;
		this.properties=[]; // array of [attribute,value] arrays (Not object to keep order)
			for (var i=6; i<data.length; i++) {this.properties.push(data[i]);}
		this.mappedValue=null;
		this.color=null;
	}

	PP.addPeptide=function(pepData) {
		peptideSet.push(new peptide(pepData));
		var lastPeptide=peptideSet[peptideSet.length-1];
		if (minValueY==null) { // 1st peptide
			minValueY=lastPeptide.valueY;
			maxValueY=lastPeptide.valueY;
		}
		else {
			minValueY=Math.min(minValueY,lastPeptide.valueY);
			maxValueY=Math.max(maxValueY,lastPeptide.valueY);
			lastAAcovered=Math.max(lastAAcovered,lastPeptide.startPos[lastPeptide.startPos.length-1]+lastPeptide.length-1);
		}
		lastPeptide.color=colorList[0]; // default
		if (peptideColor.map) {
			var mapValue;
			if (peptideColor.map==valueLabel) { // matches peptide primary value
				mapValue=convertValue(lastPeptide.value);
			}
			else if (lastPeptide[peptideColor.map.toLowerCase()] !== undefined) { // map another native peptide attribute
				mapValue=lastPeptide[peptideColor.map.toLowerCase()];
			}
			else { // Custom property
				if (mappedPropertyIdx===null) { // 1st time searched
					for (var i=0; i<peptideProperties.length; i++) {
						if (peptideProperties[i]==peptideColor.map) {
							mappedPropertyIdx=i;
							break;
						}
					}
				}
				mapValue=lastPeptide.properties[mappedPropertyIdx];
			}
			if (mapValue === undefined) {mappingError=true;}
			else {
				if (peptideColor.type.match('alpha|disc')) {
					if (!peptideColor.list) {mappedValues[mapValue]=0;} // record all possible values (no custom ordered list provided)
				}
				else if (noColorRange) { // Continuous data with no fixed range => record Min/Max as range
					if (!peptideColor.range) {
						peptideColor.range=[mapValue,mapValue];
					}
					else {
						peptideColor.range[0]=Math.min(peptideColor.range[0],mapValue);
						peptideColor.range[1]=Math.max(peptideColor.range[1],mapValue);
					}
				}
			}
			lastPeptide.mappedValue=mapValue;
		}
	}

    /********************* Chart drawing (Raphael) *********************/
	var paper; // Raphael paper
    PP.draw=function() {

		/* Canvas */
		paper=Raphael(div,canvasW,canvasH);
		/* Back panel */
		var bgPanel=paper.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		paper.rect(chartX,chartY,chartW,chartH,0).attr({fill:'#F3F3F3',stroke:'none'}); //,cursor:"crosshair"
		/* Protein axis*/
		//paper.rect(chartX,chartY+chartH,chartW,5).attr({fill:'#999999',stroke:'#999999'});
		//paper.path('M'+chartX+' '+(chartY+chartH+3)+' l'+chartW+' 0').attr({stroke:'#999999','stroke-width':6});
		/* Y axis */
		paper.path('M'+chartX0+' '+chartY+' l0 '+chartH);
		/* Axes title */
		paper.text(chartX+chartW/2,chartY+chartH+30,'Positions (aa)').attr({'font-weight':'bold','font-size':14});
		var tx=Math.min(10,chartX0-50);
		var ty=chartY+chartH/2;
		//var typeStrg=(this.type=='mean')?
		paper.text(tx,ty,chartYLabel).attr({'font-weight':'bold','font-size':14}).rotate(-90,tx,ty);

		/******* Chart plotting with default range ******/
		var minVx=1;
		var maxVx=(protein.length >= lastAAcovered)? protein.length : lastAAcovered + 20;
		var optMinVx=0;
		var tickSizeVx=(maxVx <= 100)? 10 : (maxVx <= 500)? 50 : (maxVx <= 2500)? 250 : 1000;
		var pix2valScaleX=(maxVx-minVx)/chartW;
		var optRangeY=getChartScaleRange(true,minValueY,maxValueY,chartH);
		var minVy=optRangeY[0];
		var optMinVy=optRangeY[1];
		var maxVy=optRangeY[2];
		var tickSizeVy=optRangeY[3];
		var pix2valScaleY=(maxVy-minVy)/chartH;

		/* Protein axis (drawn later in case bad prot protein length */
		if (protein.length >= lastAAcovered) { // normal case
			paper.rect(chartX0,chartY+chartH,chartW,5).attr({fill:'#999999',stroke:'#999999'});
		}
		else { // bad protein length
			var protEndPos=Math.round(protein.length/pix2valScaleX);
			paper.rect(chartX0,chartY+chartH,protEndPos,5).attr({fill:'#999999',stroke:'#999999'});
			paper.rect(chartX0+protEndPos,chartY+chartH,chartW-protEndPos,5).attr({fill:'#E0E0E0',stroke:'#E0E0E0'});
		}

		/***** Protein value=1 line *****/
		if (reference.label) {
			var pathPos=chartY0+chartH-Math.round((reference.valueY-minVy)/pix2valScaleY);
			var ref=paper.path('M'+chartX+','+pathPos+' l'+chartW+',0').attr({stroke:'#000','stroke-width':1}); //,'stroke-dasharray':'--'
			paper.path('M'+(chartX+chartW+2)+','+pathPos+' l5,-4 l0,8Z').attr({stroke:'#000',fill:'#000'});
			paper.text(chartX+chartW+9,pathPos,reference.value).attr({'font-size':11,'text-anchor':'start','font-weight':'bold',fill:'#000'});
			//Sensor
			paper.path('M'+chartX+','+pathPos+' l'+chartW+',0').attr({'stroke-width':5,'opacity':0})
			.data({posY:pathPos+0.5}) //,popupText:reference.label
			.hover(
				function(){ref.attr({'stroke-width':3}); protein.popup=drawLabel(paper,chartX,this.data('posY'),0,reference.label,'#000');},
				function(){ref.attr({'stroke-width':1}); protein.popup.remove(); protein.popup=null;}
			);
		}

		/***** Protein value & relative thresholds *****/
		if (protein.valueY != null) {
			var thList=[protein.valueY];
			if (protein.thresholdsY) {
				for (var i=0; i<protein.thresholdsY.length; i++) {thList.push(protein.thresholdsY[i]);}
			}
			var rootText=valueLabel+' '+protein.name;
			for (var i=0; i<thList.length; i++) { //for (var i=-1; i<=1; i++) {
				var pathPos=chartY0+chartH-Math.round((thList[i]-minVy)/pix2valScaleY);
				var popupText,strokeW,dashStrg;
				if (i==0) { // protein.valueY
					popupText=rootText + ': '+convertValueDisplayed(protein.value);
					strokeW=3;
					dashStrg='';
				}
				else {
					popupText=(protein.relativeThresholds[i-1] < 1)? rootText+' x 1/'+(1/protein.relativeThresholds[i-1]) : rootText+' x '+protein.relativeThresholds[i-1];
					strokeW=1;
					dashStrg='- ';
				}
				var protLine=paper.path('M'+chartX+','+pathPos+' l'+chartW+',0').attr({stroke:'#FF0000','stroke-width':strokeW,'stroke-dasharray':dashStrg});
				if (i==0) {
					paper.path('M'+(chartX+chartW+2)+','+pathPos+' l5,-4 l0,8Z').attr({stroke:'#FF0000',fill:'#FF0000'});
					paper.text(chartX+chartW+9,pathPos,convertValueDisplayed(protein.value)).attr({'font-size':11,'text-anchor':'start','font-weight':'bold',fill:'#FF0000'});
				}
				//Sensor
				paper.path('M'+chartX+','+pathPos+' l'+chartW+',0').attr({'stroke-width':5,'opacity':0})
				.data({posY:pathPos+0.5,popupText:popupText,strokeWidth:strokeW,line:protLine})
				.hover(
					function(){this.data('line').attr({'stroke-width':3}); protein.popup=drawLabel(paper,chartX,this.data('posY'),0,this.data('popupText'),'#FF0000');},
					function(){this.data('line').attr({'stroke-width':this.data('strokeWidth')}); protein.popup.remove(); protein.popup=null;}
				);
			}
		}

		/***** Axes ticks *****/
		var posX,posY;
		// X ticks
		posY=chartY+chartH+5; // 11 px below chart
		//var tickSizeVx=optRangeVx[2];
		var fixedX=(tickSizeVx < 1)? (Math.round(1/tickSizeVx)+'').length : (Math.round(tickSizeVx)==tickSizeVx)? 0 : 1;
		posX=chartX;
		paper.path('M'+posX+' '+posY+' l0 5');
		paper.text(posX,posY+10,1).attr({'font-size':10});
		var tickVx=optMinVx;
		while (tickVx <= maxVx) {
			posX=chartX0+Math.round((tickVx-minVx)/pix2valScaleX);
			if (posX >= chartX) {
				paper.path('M'+posX+' '+posY+' l0 5');
				paper.text(posX,posY+10,tickVx.toFixed(fixedX)).attr({'font-size':10});
			}
			tickVx+=tickSizeVx;
		}
		// Y ticks
		posX=chartX0-6;
		//var tickSizeVy=optRangeVy[2];
		var fixedY=(tickSizeVy < 1)? (Math.round(1/tickSizeVy)+'').length : (Math.round(tickSizeVy)==tickSizeVy)? 0 : 1;
		var tickVy=optMinVy;
		while (tickVy <= maxVy) {
			posY=chartY0+chartH-Math.round((tickVy-minVy)/pix2valScaleY);
			if (posY <= chartY0+chartH) {
				paper.path('M'+posX+' '+posY+' l5 0');
				paper.text(posX-1,posY,tickVy.toFixed(fixedY)).attr({'font-size':10,'text-anchor':'end'});
				//.hover(
				//	function(){tickPopup=drawTickPopup(this);}, // drawLabel(paper,this.attr('x'),this.attr('y'),0,'Hello','000')
				//	function(){tickPopup.remove();}
				//);
			}
			tickVy+=tickSizeVy;
		}

		/***** Plotting peptides *****/
		/* Peptide Color */
		if (peptideColor.map) {
			if (mappingError) {
				paper.text(chartX0+25,chartY0+(chartH/2),'Error: Property "'+peptideColor.map+'" not found in data!').attr({'font-size':18,fill:'#FF0000','text-anchor':'start'});
				return;
			}
			/* Plotting color palette */
			if (peptideColor.range) {
				var txt=(peptideColor.map==valueLabel)? chartYLabel : peptideColor.map;
				var t=paper.text(chartX0,12,txt+':').attr({'font-size':11,'text-anchor':'start','font-weight':'bold',fill:'#000'});
				var tBox=t.getBBox();
				var startX=chartX0 + tBox.width + 5, widthX=Math.min(200,chartW-50);
				var t2=paper.text(startX-3,12,peptideColor.range[0].toPrecision(3)).attr({'font-size':11,'text-anchor':'start',fill:'#000'});
				var t2Box=t2.getBBox();
				startX+=t2Box.width;
				var gradient=(peptideColor.order && peptideColor.order.match('desc'))? '0-rgb(0%,0%,100%)-rgb(70%,70%,0%)' : '0-rgb(70%,70%,0%)-rgb(0%,0%,100%)';
				paper.rect(startX,8,widthX,10).attr({stroke:'none',fill:gradient});
				paper.text(startX+widthX+3,12,peptideColor.range[1].toPrecision(3)).attr({'font-size':11,'text-anchor':'start',fill:'#000'});
			}
			/* Computing peptide color */
			if (peptideColor.type.match('alpha|disc')) {
				var sortedKeys;
				if (peptideColor.list) { // custom ordered list of all possible discrete values (some may not actually be used)
					sortedKeys=peptideColor.list;
				}
				else {
					var keys = Object.keys(mappedValues)
					if (peptideColor.type.match('^num')) { // numerical values
						if (peptideColor.order && peptideColor.order.match('desc')) { // descending sort
							sortedKeys=keys.sort(function(a,b){return b-a});
						}
						else { // ascending sort
							sortedKeys=keys.sort(function(a,b){return a-b});
						}
					}
					else { // alpha-numerical values
						sortedKeys=keys.sort();
						if (peptideColor.order && peptideColor.order.match('desc')) {sortedKeys.reverse();}
					}
				}
				/* Record order in mappedValues object */
				for (var i=0; i<sortedKeys.length; i++) {
					mappedValues[sortedKeys[i]]=i;
				}
				/* Apply discrete color mapping to peptides */
				for (var i=0; i<peptideSet.length; i++) {
					var pepValueIdx=mappedValues[peptideSet[i].mappedValue] % colorList.length;
					peptideSet[i].color=colorList[pepValueIdx];
				}
				/* Draw list of discrete values color & name */
				var x=(sortedKeys.length >= 10)? 5: chartX0;
				var yb=8,yt=12;
				for (var i=0; i<sortedKeys.length; i++) {
					if (i==10) {
						var t=paper.text(x,yt,'...').attr({'font-size':11,'text-anchor':'start','font-weight':'bold',fill:'#000'});
						x+=t.getBBox().width + 10;
						break;
					}
					var idx=i % colorList.length;
					paper.rect(x,yb,10,10).attr({stroke:'none',fill:colorList[idx]});
					x+=12;
					var t=paper.text(x,yt,sortedKeys[i]).attr({'font-size':11,'text-anchor':'start','font-weight':'bold',fill:'#000'});
					x+=t.getBBox().width + 10;
				}
				if (x > canvasW) { // increase paper width to fit all legends
					canvasW=x;
					paper.setSize(canvasW,canvasH);
					bgPanel.attr({width:canvasW});
				}
			}
			else { // continuous values
				var delta=peptideColor.range[1]-peptideColor.range[0];
				for (var i=0; i<peptideSet.length; i++) {
					var relValueDown,relValueUp;
					if (peptideColor.order && peptideColor.order.match('desc')) {
						relValueDown=Math.round(70*(Math.max(0,Math.min(delta,peptideSet[i].mappedValue-peptideColor.range[0])))/delta);
						relValueUp=Math.round(100*(Math.max(0,Math.min(delta,peptideColor.range[1]-peptideSet[i].mappedValue)))/delta);
					}
					else { // asc
						relValueDown=Math.round(70*(Math.max(0,Math.min(delta,peptideColor.range[1]-peptideSet[i].mappedValue)))/delta);
						relValueUp=Math.round(100*(Math.max(0,Math.min(delta,peptideSet[i].mappedValue-peptideColor.range[0])))/delta);
					}
					peptideSet[i].color='rgb('+relValueDown+'%,'+relValueDown+'%,'+relValueUp+'%)'; // 70%,70%,0% <-> 0%,0%,100% (dark yellow <-> blue)
				}
			}
		}
		/* Peptide coverage */
		for (var i=0; i<peptideSet.length; i++) {
			var pSizeX=Math.round(peptideSet[i].length/pix2valScaleX)+2; // +2 to compensate for invisible stroke width
			for (var j=0; j<peptideSet[i].startPos.length; j++) { // multiple occurences
				var pX=chartX0+Math.round((peptideSet[i].startPos[j]-minVx)/pix2valScaleX)-1; // -1 for invisible stroke width
				paper.rect(pX,chartY+chartH,pSizeX,5).attr({stroke:'none',fill:'#505050'});
			}
		}
		/* Peptides & protein projections */
		for (var i=0; i<peptideSet.length; i++) {
			var pX=Math.round((peptideSet[i].startPos[0]-minVx)/pix2valScaleX)+chartX0-1.5; // -1 for invisible stroke width
			var pY=chartY0+chartH-Math.round((peptideSet[i].valueY-minVy)/pix2valScaleY)-2.5; // -2 for invisible stroke width
			var pSizeX=Math.round(peptideSet[i].length/pix2valScaleX)+2; // +2 to compensate for invisible stroke width
			var p;
			if (peptideSet[i].excluded) { // excluded peptide
				var pathStrg=(peptideSet[i].varMod)? 'M'+(pX+1)+' '+(pY-1)+'l'+(pSizeX-2)+' 0 l0 3 l'+(2-pSizeX+3)+' 0 l-3 3 Z' : 'M'+(pX+1)+' '+(pY-1)+'l'+(pSizeX-2)+' 0 l0 3 l'+(2-pSizeX)+' 0 Z';
				p=paper.path(pathStrg).attr({stroke:peptideSet[i].color,'stroke-opacity':1,fill:peptideSet[i].color,'fill-opacity':0});
			}
			else {
				var pathStrg=(peptideSet[i].varMod)? 'M'+pX+' '+pY+'l'+pSizeX+' 0 l0 5 l'+(-pSizeX+4)+' 0 l-4 4 Z' : 'M'+pX+' '+pY+'l'+pSizeX+' 0 l0 5 l'+(-pSizeX)+' 0 Z';
				p=paper.path(pathStrg).attr({stroke:peptideSet[i].color,'stroke-opacity':0,fill:peptideSet[i].color,'fill-opacity':0.75});
			}
			p.data({peptide:peptideSet[i],showLabel:null}) //x:pX,y:pY,
			.hover(function(){emphasizePeptide(this,'on')},function(){emphasizePeptide(this,'off')})
			.click(function(){selectPeptide(this,'auto')});
			var pSet=paper.set();
			for (var j=0; j<peptideSet[i].startPos.length; j++) { // multiple occurences
				var pXj=Math.round((peptideSet[i].startPos[j]-minVx)/pix2valScaleX)+chartX0-1;
				pSet.push(paper.rect(pXj,chartY+chartH,pSizeX,5));
			}
			p.data({proteinProj:pSet.attr({stroke:'none',fill:'#FF0000'}).hide()});
		}
	}

    /******************* Peptide selection & highlight ***********************/
	function getPeptideInfo(dp) {
		var infoStrg=dp.sequence;
		if (dp.varMod) infoStrg+=' + '+dp.varMod;
		if (dp.charge) infoStrg+='\nCharge: '+dp.charge+'+';
		infoStrg+='\nPos: ';
		for (var i=0; i<dp.startPos.length; i++) {
			if (i>0) infoStrg+=',';
			infoStrg+=dp.startPos[i]+'-'+(dp.startPos[i]*1+dp.length-1);
		}
		infoStrg+='\n'+valueLabel+': '+convertValueDisplayed(dp.value);
		if (dp.valueStrg) {infoStrg+=' ['+dp.valueStrg+']';}
		if (dp.excluded) {infoStrg+='\n***Excluded***';}
		for (var i=0; i<peptideProperties.length; i++) {
			if (peptideProperties[i]=='id') continue; // skip id if any
			if (dp.properties[i] != null) infoStrg+='\n'+peptideProperties[i]+': '+dp.properties[i];
		}
		return infoStrg;

	}
	function selectPeptide(p,action) {
		if (action=='auto') { // invert selection
		    if (p.data('showLabel')) { // already showing label => set OFF
			    p.data('showLabel',null);
			    emphasizePeptide(p,'off');
		    }
		    else { // set ON
			    p.data('showLabel',1);
			    emphasizePeptide(p,'on');
		    }
	    }
	    else if (action=='on') { // set ON
		    if (!p.data('showLabel')) {
			    p.data('showLabel',1);
			    emphasizePeptide(p,'on');
		    }
	    }
	    else { // OFF
		    p.data('showLabel',null);
		    emphasizePeptide(p,'off');
	    }
	}
	function emphasizePeptide(p,action) {
		var pColor;
		var dp=p.data('peptide');
		if (action=='on') {
			pColor='#FF0000';
			if (p.data('labelObj')) {p.data('labelObj').show();}
			else {displayPeptideLabel(p);}
			p.data('proteinProj').show();
		}
		else { // off
			if (p.data('showLabel')) { // needed when hovering out of point
				pColor='#FF0000';
			}
			else {
				pColor=dp.color;
				if (p.data('labelObj')) {
					p.data('labelObj').remove();
					p.removeData('labelObj');
				}
				p.data('proteinProj').hide();
			}
		}
		p.attr({stroke:pColor,fill:pColor});
    }
    function displayPeptideLabel(p) {
	    var dp=p.data('peptide');
		var pBox=p.getBBox();
	    var text=getPeptideInfo(dp);
	    //var tSet=drawPeptideLabel(p.attr('x'),p.attr('y'),p.attr('width'),text,lColor);
	    //tSet[0].data({x:p.attr('x'),y:p.attr('y')});
		var tSet=drawPeptideLabel(pBox.x-0.5,pBox.y,pBox.width,text,dp.color);
		if (plotData.peptideOnClick) {
			tSet.click(function(){plotData.peptideOnClick(dp.properties[pepIdPropertyIndex],dp)})
			.attr({cursor:'pointer'});
		}
	    tSet[0].data({x:pBox.x,y:pBox.y});
//tSet.translate(0,50);
//console.log('SET: x='+tSet[0].data('x')+', y='+tSet[0].data('y'));
	    p.data('labelObj',tSet);
    }
    function drawPeptideLabel(x,y,pSize,text,lColor) {
//console.log(lColor);
	    var shift=15;
	    var d=0;
	    var t=paper.text(x+pSize+shift,y-shift,text).attr({'font-size':10,'text-anchor':'start',fill:lColor}); //,'font-weight':'bold','fill-opacity':0.6
	    var tBox=t.getBBox();
	    var tx=tBox.x;
	    var ty=tBox.y;
	    var tw=tBox.width;
	    var th=tBox.height;
	    //var tb; //Bubble around text
	    //var dy=Math.round((th-6)/2);
	    var th05=Math.round(th/2);
	    var sx,ex;
	    if (tx+tw > paper.width) {
		    //t.attr({x:paper.width-tBox.width});
		    tx=tx-pSize-2*(shift)-tw; //tBox.width;
		    t.attr({x:tx});
		    //tb=paper.path('M'+(tx-2)+' '+(ty-1)+' l'+(tw+4)+' 0 l0 '+dy+' l10 3 l-10 3 l0 '+dy+' l'+(-tw-4)+' 0 Z');
		    sx=tx+tw+3;
		    ex=x;
	    }
	    else {
		    //tb=paper.path('M'+(tx-2)+' '+(ty-1)+' l'+(tw+4)+' 0 l0 '+(th+2)+' l'+(-tw-4)+' 0 l0 '+(-dy)+' l-10 -3 l10 -3 Z');
		    sx=tx-1;
		    ex=x+pSize;
	    }
	    if (ty < 0) {
			t.attr({y:2+th05});
			ty=t.getBBox().y;
	    }
	    var sy=ty+th05;
	    //var ex=x+pSize; //Math.round(d*(sx-x)/(d+shift));
	    var ey=y+Math.round(d*(sy-y)/(d+shift));
	    var tl=paper.path('M'+sx+' '+sy+' L'+ex+' '+ey).attr({stroke:lColor});
	    var tb=paper.rect(tx-2,ty-1,tw+4,th+2,4).attr({stroke:lColor,fill:'#fff','fill-opacity':0.7});
	    t.toFront();

	    return paper.set(tb,tl,t);
    }
/*
	function drawTickPopup(v) {
		var log2val=v.attr('text');
		var val=Math.exp(log2val*0.693147181);
		if (val >= 1) {val=val.toFixed(2);}
		else {val='1/'+(1/val).toFixed(2);}
		var t=paper.text(v.attr('x'),v.attr('y')-17,'x'+val).attr({'font-size':12,'text-anchor':'end'});
		var tBox=t.getBBox();
	    var tx=tBox.x;
	    var ty=tBox.y;
	    var tw=tBox.width;
	    var th=tBox.height;
		var tb=paper.rect(tx-2,ty-1.5,tw+4,th+2,0).attr({stroke:'#000',fill:'#fff'});
		t.toFront();
		return paper.set(tb,t);
	}
*/
} // end of PP

/*
####>Revision history<####
# 1.1.4 Handles custom ordered list of discrete values in peptideColor.list=[array] (PP 29/03/17)
# 1.1.3 Handles bad/missing protein length (PP 16/02/17)
# 1.1.2 Now uses getChartScaleRange from chartLibrary2.js (PP 27/05/16)
# 1.1.1 Add peptideOnclick option on popup label. Requires peptide externalID. Now requires chartLibrary2.js for line popup label (PP 18/05/16)
# 1.1.0 Major update to allow visualization any kind peptide data (PP 14/04/16)
# 1.0.5 Bug fix in display of protein ratio < 1 (PP 03/09/15)
# 1.0.4 Draws protein ratio or protein mean (PP 16/07/15)
# 1.0.3 Add reference FC=1 option & minor display improvements (PP 11/12/14)
# 1.0.2 Does not display test/ref for FC in peptide popup (PP 03/10/14)
# 1.0.1 GPL license (PP 23/09/13)
# 1.0.0 Production version (PP 01/03/12)
*/
