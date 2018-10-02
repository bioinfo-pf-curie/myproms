/*
################################################################################
# spectrumPlot.js    1.0.2                                                     #
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
/*
function ionFragment(ionString) { // mz,intensity,label '627.32:0.03' '459.30:0.20=y4'
	var data=ionString.split(/[:=]/);
	this.mz=data[0];
	this.position=[0,0]; // position on chart in px
	this.intensity=data[1];
    this.label=data[2] || null; // possibly null
    this.aspect=[]; // Array of svg elements ID [ion frgmt,label/tip,connec. line]
}
*/

/*****************************************************************
                  Generic Plot object
******************************************************************/
function spectrumPlot(spData) {
//console.time('LOAD');
    var SP=this;
	//this.jsName=spData.name;
    this.divID=spData.div;
    this.div=document.getElementById(this.divID);
	this.ready=false;

    /******** Chart variables & objects ********/
	var verticalIonLabel=(spData.isVert)? true : false;
	this.isReference=(spData.isRef)? true : false;
	var ionColor=(this.isReference)? '#00F' : '#F00';
	this.ionOnclick=spData.ionOnclick;
	this.ionOnHover=spData.ionOnHover;
	this.onViewChange=spData.onViewChange;
	this.unifyView=false;
	var firstPlot=true;
	var ionPopup; // svg popup

	this.zoomable=true;
	this.autoZoom=true; // immediate zoom after area selection
	this.hideZoom=true;
	this.noZoomY=true;

	this.chartSettings={reference:{},current:{}};
	this.chartMarks=[];
	this.plotArea;
	this.dragArea;
	this.dragContext='area';
	//this.zoomText;
	this.panProcess={x:0,y:0,dx:0,dy:0}; //new Object();
	this.highlightArea;

	/***** Chart axes parameters *****/
	this.axisXtext='m/z';
	this.axisYtext=''; //'Intensity %';
this.noAxisY=true;
	this.noTicksY=true;
	//this.forceXto0=1; // 0,1,2 (not true,false!)
//this.keepXabove0=true;
	this.forceYto0=1; // 0,1,2 (not true,false!)

	/******** Chart geometry ********/
    var BORDER_SPACE=10;
    var XAXIS_SPACE=40;
    var YAXIS_SPACE=0; //50;
	var LABEL_SPACE=0; //(verticalIonLabel)? 20 : 10;
    var chartW=(spData.width)? spData.width : 800;
    var chartH=(spData.height)? spData.height : 200;
    var chartX=BORDER_SPACE+YAXIS_SPACE+1;
    var chartY=BORDER_SPACE+LABEL_SPACE+1;
	var canvasW=chartX+chartW+BORDER_SPACE-1;
    var canvasH=chartY+chartH+XAXIS_SPACE+BORDER_SPACE-1;
	var baseLine=chartY+chartH;
	var scrollLeftPath='M16,30.534c8.027,0,14.534-6.507,14.534-14.534c0-8.027-6.507-14.534-14.534-14.534C7.973,1.466,1.466,7.973,1.466,16C1.466,24.027,7.973,30.534,16,30.534zM18.335,6.276l3.536,3.538l-6.187,6.187l6.187,6.187l-3.536,3.537l-9.723-9.724L18.335,6.276z';
	var scrollRightPath='M16,1.466C7.973,1.466,1.466,7.973,1.466,16c0,8.027,6.507,14.534,14.534,14.534c8.027,0,14.534-6.507,14.534-14.534C30.534,7.973,24.027,1.466,16,1.466zM13.665,25.725l-3.536-3.539l6.187-6.187l-6.187-6.187l3.536-3.536l9.724,9.723L13.665,25.725z';
	var zoomMinusPath='M22.646,19.307c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127l3.535-3.537L22.646,19.307zM13.688,20.369c-3.582-0.008-6.478-2.904-6.484-6.484c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486C20.165,17.465,17.267,20.361,13.688,20.369zM8.854,11.884v4.001l9.665-0.001v-3.999L8.854,11.884z';
	var zoomOutPath='M29.772,26.433l-7.126-7.126c0.96-1.583,1.523-3.435,1.524-5.421C24.169,8.093,19.478,3.401,13.688,3.399C7.897,3.401,3.204,8.093,3.204,13.885c0,5.789,4.693,10.481,10.484,10.481c1.987,0,3.839-0.563,5.422-1.523l7.128,7.127L29.772,26.433zM7.203,13.885c0.006-3.582,2.903-6.478,6.484-6.486c3.579,0.008,6.478,2.904,6.484,6.486c-0.007,3.58-2.905,6.476-6.484,6.484C10.106,20.361,7.209,17.465,7.203,13.885z';

    /********************* Data import *********************/
//this.dataSets=[];
	var tolerance=spData.tolerance || 0; // tolerance on delta theo vs obs mz
	this.ionSeries=[];
	this.minValueX=2000;
	this.minValueY=this.maxValueX=this.maxValueY=0;

    /***** Ion series *****/
	var ionSet=spData.ionSeries.split(';');
	for (var i=0; i < ionSet.length; i++) {
		if (!ionSet[i]) continue; // tailing ';'
		//var ion=new ionFragment(ionSet[i]);
		var data=ionSet[i].split(/[:=]/);
		var ion={mz:data[0],
				position:[0,0], // position on chart in px
				intensity:data[1],
				label:data[2] || null, // possibly null
				aspect:[] // [fragment,label,link] link is optional
				};
		this.ionSeries.push(ion);
		this.minValueX=Math.min(this.minValueX,ion.mz);
		this.maxValueX=Math.max(this.maxValueX,ion.mz);
		this.maxValueY=Math.max(this.maxValueY,ion.intensity);
    }
	/*Converting intensities to % */
	this.maxValueY/=100;
	for (var i=0; i < this.ionSeries.length; i++) {
		this.ionSeries[i].intensity/=this.maxValueY;
	}
	this.maxValueY=100;
//console.timeEnd('LOAD');

	/********************* Plot display (Raphael) *********************/
    this.draw=function() {
		/* Canvas */
		this.canvas=Raphael(this.divID,canvasW,canvasH);
//this.canvas.renderfix();
		/* Back panel */
		this.backPanel=this.canvas.rect(0,0,canvasW,canvasH,0).attr({fill:'#fff',stroke:'#000'});
		/* Chart area */
		this.plotArea=this.canvas.rect(chartX,chartY,chartW,chartH,0)
		.attr({fill:'#F3F3F3',stroke:'none'}) //,cursor:"crosshair"
		/*.hover(
			function() {
				if (SP.chartSettings.current.X.minRange > SP.chartSettings.reference.X.minRange) {SP.scrollLeftIcon.toFront().show();}
				if (SP.chartSettings.current.X.maxRange < SP.chartSettings.reference.X.maxRange) {SP.scrollRightIcon.toFront().show();}
			},
			function() {SP.scrollLeftIcon.hide(); SP.scrollRightIcon.hide();}
		)*/
		.drag(function(dx,dy,x,y,e) {extendDragging(SP,dx,dy,x,y,e);},
			  function(x,y,e) {setDragging(SP,x,y,e);},
			  function() {endDragging(SP);}
			  )
		.dblclick(function(){zoomOut(SP)});
		if (LABEL_SPACE) {
			this.canvas.rect(chartX,BORDER_SPACE+1,chartW,LABEL_SPACE,0)
			.attr({fill:'#F3F3F3',stroke:'none'});
		}

		/*** highlightArea ***/
		this.highlightArea=this.canvas.set(
			this.canvas.rect(chartX,BORDER_SPACE+1,3,chartH+LABEL_SPACE,0).attr({fill:'#0FF',stroke:'none','fill-opacity':0.5}),
			this.canvas.path('M0 '+(BORDER_SPACE+1)+' l0 '+(chartH+LABEL_SPACE)).attr({stroke:'#099','stroke-dasharray':'-'})
			).hide();

		/***** Display chart(s) *****/
		initializeChart(SP);

		this.scrollLeftIcon=this.canvas.path(scrollLeftPath).attr({fill:'#000',stroke:'none','fill-opacity':0.3}).transform('T'+(chartX+5)+','+(chartY+(chartH/2)-14)).hide()
							.hover(function(){this.show(); this.attr({fill:'#D00','fill-opacity':1});},function(){this.attr({fill:'#000','fill-opacity':0.3});})
							.mousedown(function(){SP.scrollChart(1)});
		this.scrollRightIcon=this.canvas.path(scrollRightPath).attr({fill:'#000',stroke:'none','fill-opacity':0.3}).transform('T'+(chartX+chartW-37)+','+(chartY+(chartH/2)-14)).hide()
							.hover(function(){this.show(); this.attr({fill:'#D00','fill-opacity':1});},function(){this.attr({fill:'#000','fill-opacity':0.3});})
							.click(function(){SP.scrollChart(-1)});

		this.zoomMinusIcon=this.canvas.path(zoomMinusPath).attr({fill:'#000',stroke:'none','fill-opacity':0.3}).transform('T'+chartX+','+(chartY+chartH+10)).hide();
		var zoomMinArea=this.canvas.rect(chartX+3,chartY+chartH+12,25,28).attr({stroke:'none',fill:'#000','fill-opacity':0})
							.hover(function(){SP.zoomMinusIcon.attr({fill:'#D00','fill-opacity':1});},function(){SP.zoomMinusIcon.attr({fill:'#000','fill-opacity':0.3});})
							.click(function(){zoomOut(SP)});

		this.zoomOutIcon=this.canvas.path(zoomOutPath).attr({fill:'#000',stroke:'none','fill-opacity':0.3}).transform('T'+(chartX+30)+','+(chartY+chartH+10)).hide();
		var zoomOutArea=this.canvas.rect(chartX+33,chartY+chartH+12,25,28).attr({stroke:'none',fill:'#000','fill-opacity':0})
							.hover(function(){SP.zoomOutIcon.attr({fill:'#D00','fill-opacity':1});},function(){SP.zoomOutIcon.attr({fill:'#000','fill-opacity':0.3});})
							.click(function(){zoomOut(SP,true)});

		this.ready=true;
	}

	this.plotChart=function() {
//console.time('DRAW_INIT');
		var settings=this.chartSettings.current;
		if (this.unifyView==true) {this.onViewChange(this);}
		else if (this.onViewChange) {this.unifyView=true;} // for future view changes

		/***** Ions & labels *****/
		var ionScale=(LABEL_SPACE)? 1 : (verticalIonLabel)? 0.9 : 0.95;
		var labelListIdx=[];
		for (var i=0; i<this.ionSeries.length; i++) {
			var ion=this.ionSeries[i];
			if (ion.mz >= settings.X.minRange && ion.mz <= settings.X.maxRange) { // within current range
				var iX=chartX+Math.round((ion.mz-settings.X.minRange)/settings.X.pix2valRatio)+0.5;
				var iI;
				if (firstPlot==true) { // draw fragments only once!
					iI=Math.round(ionScale*ion.intensity/settings.Y.pix2valRatio);
					ion.position[1]=iI;
					ion.aspect[0]=this.canvas.path('M'+iX+' '+baseLine+' l0 -'+iI).attr({stroke:ionColor}).data({ionIdx:i})
									.hover(function(){_displayMZ(this,'on')},function(){_displayMZ(this,'off')});
					if (ion.label) {
						if (ion.label==' ') { // matched but no displayed label
							ion.aspect[1]=this.canvas.path('M'+iX+' '+(baseLine-iI-2)+' l0 -7').data({ionIdx:i}) //.attr({stroke:'#F00'});
									.hover(function(){_displayMZ(this,'on')},function(){_displayMZ(this,'off')});
						}
						else { // displayed label
							labelListIdx.push([i,iI]);
						}
					}
					ion.position[0]=iX;
				}
				else { // translate
					var dx=iX-ion.position[0];
					for (var j=0; j<ion.aspect.length; j++) {
						//ion.aspect[j].translate(dx,0).show();
						var transformStrg='t'+dx+',0';
						if (ion.aspect[j].type=='text' && verticalIonLabel) {transformStrg+='r-90,'+ion.aspect[j].attr('x')+','+ion.aspect[j].attr('y');}
						ion.aspect[j].transform(transformStrg).show();
					}
					iI=ion.position[1];
				}
				//ion.position[0]=iX;
			}
		}
//console.timeEnd('DRAW_INIT');
//console.time('LABELS #1');
		var boxList=[]; // list of text boxes for overlap check
		var labelShift=30;
		for (var i=0; i<labelListIdx.length; i++) { // ****empty if ! firstPlot****
			var ion=this.ionSeries[labelListIdx[i][0]];
			var label=ion.label;
			var iX=ion.position[0];
			var iX2=iX;
			var iI=labelListIdx[i][1];
			var lY=(verticalIonLabel)? baseLine-iI-3 : baseLine-iI-7;
			var movedLabel=false;
			if (lY > baseLine-labelShift) {
				lY=baseLine-labelShift;
				movedLabel=true;
			}
			var l=this.canvas.text(iX,lY,label).attr({'font-size':12,'text-anchor':'start'}); //,fill:'#00F'
			if (verticalIonLabel) {
				//l.rotate(-90,iX,lY);
				l.transform('r-90,'+iX+','+lY);
			}
			l.data({ionIdx:labelListIdx[i][0]})
			.hover(function(){_displayMZ(this,'on')},function(){_displayMZ(this,'off')});

			// Checking for label overlap
			var lBox=l.getBBox();
			var overlap=true; //false; //
			while (overlap) {
				overlap=false;
				for (var b=0; b<boxList.length; b++) {
					var lastLoop=false;
					if (Raphael.isBBoxIntersect(lBox,boxList[b])) { // while
						overlap=true; // intersect event
						lY=boxList[b].y-lBox.height-5; // make sure no more overlap
						if (lY < 30) {
							lY=30;
							// Switch to horizontal shift
							/*
							if (lBox.x >= boxList[b].x) {
								iX2=boxList[b].x2+5;
							}
							else {iX2=boxList[b].x-lBox.width-5;}
							l.attr({x:iX2});
							*/
							lastLoop=true;
						}
						l.attr({y:lY});
						if (verticalIonLabel) {l.transform('r-90,'+l.attr('x')+','+l.attr('y'));} // update rotation to match new y
						lBox=l.getBBox();
						movedLabel=true;
						if (lastLoop) {
							overlap=false;
							break;
						}
					}
				}
			}
			if (verticalIonLabel) {
				if (lBox.y < 0) { // while
					var txt=l.attr('text').substring(0,4);
					l.attr({text:txt+'>'});
					lBox=l.getBBox();
				}
			}
			else { // horizontal
				if (lBox.x+lBox.width > canvasW) { // while
					var txt=l.attr('text').substring(0,4);
					l.attr({text:txt+'>'});
					lBox=l.getBBox();
				}
			}
			boxList.push(lBox);
			ion.aspect[1]=l;
			// Adding connection line
			if (movedLabel) {
				var cl=this.canvas.path('M'+iX+' '+(baseLine-iI-2)+' L'+iX2+' '+lBox.y2).attr({'stroke-dasharray':'.',stroke:'#000'});
				ion.aspect[2]=cl;
			}
			labelShift+=10;
			if (labelShift > 100) {labelShift=30;}
		}
//console.timeEnd('LABELS #1');

		/***** Scroll/zoom icons *****/
		if (settings.X.minRange > this.chartSettings.reference.X.minRange) {
			this.scrollLeftIcon.toFront().show();
		}
		if (settings.X.maxRange < this.chartSettings.reference.X.maxRange) {
			this.scrollRightIcon.toFront().show();
		}
		if (settings.X.minRange > this.chartSettings.reference.X.minRange || settings.X.maxRange < this.chartSettings.reference.X.maxRange) {
			this.zoomMinusIcon.show();
			this.zoomOutIcon.show();
		}

		/***** highlightArea width *****/
		var hw=Math.max(3,Math.round(2*tolerance/settings.X.pix2valRatio));
		this.highlightArea[0].attr({'width':hw});

		firstPlot=false;
    }

	this.clearChart=function() {
		for (var i=0; i<this.ionSeries.length; i++) {
			for (var j=0; j<this.ionSeries[i].aspect.length; j++) { // clear label
				this.ionSeries[i].aspect[j].hide();
			}
		}
		this.scrollLeftIcon.hide();
		this.scrollRightIcon.hide();
		this.zoomMinusIcon.hide();
		this.zoomOutIcon.hide();
		viewHasChanged=true; // flag to tell a change occured in view (zoom/scroll)
	}

	this.scrollChart=function(sens)  {
		this.panProcess.dx=sens*250;
		this.panProcess.dy=0;
		panChart(this);
	}

	this.highlightMZ=function(action,mz) {
		if (mz >= this.chartSettings.current.X.minRange && mz <= this.chartSettings.current.X.maxRange) { // in range
			if (action=='on') {
				var hX=chartX+Math.round((mz-this.chartSettings.current.X.minRange)/this.chartSettings.current.X.pix2valRatio)+0.5;
				this.highlightArea[0].attr({x:hX-this.highlightArea[0].attr('width')/2});
				this.highlightArea[1].transform('T'+hX+',0');
				this.highlightArea.show();
			}
			else {
				this.highlightArea.hide();
			}
		}
		else if (mz < this.chartSettings.current.X.minRange) {
			if (action=='on') {
				this.scrollLeftIcon.attr({fill:'#0FF','fill-opacity':0.5}).show();
			}
			else {
				this.scrollLeftIcon.attr({fill:'#000','fill-opacity':0.3});
				if (this.chartSettings.current.X.minRange == this.chartSettings.reference.X.minRange) {this.scrollLeftIcon.hide();}
			}
		}
		else {
			if (action=='on') {
				this.scrollRightIcon.attr({fill:'#0FF','fill-opacity':0.5}).show();
			}
			else {
				this.scrollRightIcon.attr({fill:'#000','fill-opacity':0.3});
				if (this.chartSettings.current.X.maxRange == this.chartSettings.reference.X.maxRange) {this.scrollRightIcon.hide();}
			}
		}
	}

	function _displayMZ(el,action) {
		if (action=='on') {
			var label=SP.ionSeries[el.data('ionIdx')].label;
			var mz=SP.ionSeries[el.data('ionIdx')].mz;
			var popStrg=(label && label != ' ')? 'Ion: '+label+'\nm/z: '+mz : 'm/z: '+mz;
			var box=el.getBBox();
			var x=box.x+(box.x2-box.x)/2;
			ionPopup=drawLabel(SP.canvas,x,box.y,0,popStrg,'#000');
		}
		else {
			ionPopup.remove();
		}
	}

	this.toBlackAndWhite=function(bw) {
		if (bw) {
			this.plotArea.attr({fill:'#FFFFFF'});
			ionColor='#000';
		}
		else {
			this.plotArea.attr({fill:'#F3F3F3'});
			ionColor=(this.isReference)? '#00F' : '#F00';
		}
		for (var i=0; i<this.ionSeries.length; i++) {
			this.ionSeries[i].aspect[0].attr({stroke:ionColor});
		}
	}
} // end of SP

/*
####>Revision history<####
# 1.0.2 Added black and white / colored switch function (PP 19/06/14)
# 1.0.1 Fix vertical label zoom + label truncation when out of canvas (PP 16/04/14)
# 1.0.0 Started (PP 31/03/14)
*/
