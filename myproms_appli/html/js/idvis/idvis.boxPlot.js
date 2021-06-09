/*
################################################################################
# idvis.boxPlot.js    2.0.3                                                    #
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
                  Box object
******************************************************************/
idvis.boxPlot = function(bpData) {
	var BP=this;
	this.chartID=idvis.lib.registerChart(this);

	/* Configuration variables */
    var mainDivID=bpData.div; //document.getElementById(bpData.div);
	this.convertValue=function(axis,value) {return value;}; // for threshold line
	var horizontalBP=(bpData.horizontal)? true : false;
	this.forceTo0=bpData.forceTo0 || null;
    var canvasDivID,formDivID;
	var boxThresholdLine,editableBP;
    if (bpData.threshold && bpData.threshold.editable) { // null if not declared
		editableBP=true;
		canvasDivID=mainDivID+'_canvas';
		formDivID=mainDivID+'_form';
    }
    else {
		editableBP=false;
		canvasDivID=mainDivID;
    }
	this.getDivID=function() {return canvasDivID;};
	if (bpData.threshold) {
		let thAxis=(horizontalBP)? 'X' : 'Y';
		//boxThresholdLine=new idvis.thresholdLine(BP,thAxis,bpData.threshold.name,bpData.threshold.value,'#F00','--');
        boxThresholdLine=new idvis.feature(BP,{type:'line',axis:thAxis,label:bpData.threshold.name,value:bpData.threshold.value,color:'#F00',pattern:'- ',editable:bpData.threshold.editable});
	}

    var boxOnClick=(bpData.boxOnClick)? bpData.boxOnClick : null;
	var excludeOnClick=(bpData.excludeOnClick)? bpData.excludeOnClick : null;
	//var outlierOnClick=(bpData.outlierOnClick)? bpData.outlierOnClick : function(){};

	/*** Data variables ***/
    var toScale=(bpData.toScale)? bpData.toScale : function(value) {return value;};
	this.outlierRule=(idvis.isNumber(bpData.outlierRule))? bpData.outlierRule : 'IQR';
	var boxList=[];

	/********** Add a single box ***********/
	var maxLabelLength,minBoxPlotValue,maxBoxPlotValue;
	this.addBox=function(boxInfo) {
		let newBox=new box(this,boxList.length,boxInfo);
		boxList.push(newBox);
		if (newBox.index===0) {
			maxLabelLength=newBox.label.length;
			if (boxThresholdLine) {
				minBoxPlotValue=Math.min(boxThresholdLine.initValue,newBox.extremeValues[0]);
				maxBoxPlotValue=Math.max(boxThresholdLine.initValue,newBox.extremeValues[1]);
			}
			else {
				minBoxPlotValue=newBox.extremeValues[0];
				maxBoxPlotValue=newBox.extremeValues[1];
			}
		}
		else {
			maxLabelLength=Math.max(maxLabelLength,newBox.label.length);
			minBoxPlotValue=Math.min(minBoxPlotValue,newBox.extremeValues[0]);
			maxBoxPlotValue=Math.max(maxBoxPlotValue,newBox.extremeValues[1]);
		}
	};

    /*** Graphics variables ***/
	const topSpace=20,
		  bottomSpace=20,
		  leftSpace=20,
		  rightSpace=20; //+(lastColLabelSize*5);
	var highlightZone,
		minValue,maxValue,pix2valScale,thresholdValue,thresholdInputID,
		canvasWidth,canvasHeight,boxPlotX0,boxPlotY0,boxPlotW,boxPlotH,boxSpace,boxSize,labelFontSize;

	/********************* Box plot display (Raphael) *********************/
    this.draw=function() {
//console.log('Load='+(new Date().getTime()-t0));
		/* Computing Graphic area size */
		if (horizontalBP) {
			if (bpData.boxSize) {
				boxSize=Math.max(7,bpData.boxSize);
				boxSpace=Math.round(boxSize*1.2);
				boxPlotH=boxSpace*boxList.length;
			}
			else if (bpData.height) {
				boxPlotH=bpData.height;
				boxSpace=bowPlotH/boxList.length;
				boxSize=boxSpace/1.2;
			}
			else { // auto
				boxSpace=Math.min(150,Math.max(11,750/boxList.length));
				boxSize=Math.round(boxSpace/1.2);
				boxPlotH=boxSpace*boxList.length;
			}
			boxPlotW=(bpData.width)? bpData.width : 250;
			labelFontSize=(boxSpace > 20)? 14 : 10;
			let labelSpaceX=5+Math.round(maxLabelLength*labelFontSize/1.6);
			boxPlotX0=leftSpace+labelSpaceX+1;
			boxPlotY0=topSpace+1;
			canvasWidth=leftSpace+labelSpaceX+boxPlotW+rightSpace;
			canvasHeight=topSpace+boxPlotH+30+bottomSpace;
		}
		else { // vertical (normal case)
			if (bpData.boxSize) {
				boxSize=Math.max(7,bpData.boxSize);
				boxSpace=Math.round(boxSize*1.2);
				boxPlotW=boxSpace*boxList.length;
			}
			else if (bpData.width) {
				boxPlotW=bpData.width;
				boxSpace=boxPlotW/boxList.length;
				boxSize=boxSpace/1.2;
			}
			else { // auto
				boxSpace=Math.min(150,Math.max(11,750/boxList.length));
				boxSize=Math.round(boxSpace/1.2);
				boxPlotW=boxSpace*boxList.length;
			}
			boxPlotH=(bpData.height)? bpData.height : 250;
			labelFontSize=(boxSpace > 20)? 14 : 10;
			let labelSpace=5+Math.round(maxLabelLength*labelFontSize/1.6);
			boxPlotX0=leftSpace+50.5; //+Math.max(30,labelSpace)+1;
			boxPlotY0=topSpace+1.5;
			canvasWidth=boxPlotX0-1+boxPlotW+rightSpace;
			canvasHeight=topSpace+boxPlotH+labelSpace+bottomSpace;
		}
//console.log(canvasWidth+' * '+canvasHeight);
		/* DIVs */
		if (editableBP) {
		    //document.getElementById(mainDivID).innerHTML='<TABLE><TR><TD valign="top"><DIV id="'+formDivID+'"></DIV></TD></TR><TR><TD><DIV id="'+canvasDivID+'"></DIV></TD></TR></TABLE>';
			idvis.lib.drawChartLayout(BP,bpData.formPosition,canvasWidth,canvasDivID,formDivID);

			/* Form menu */
			//initializeForm();
			thresholdInputID='thresValueID_'+BP.chartID;
			let htlmCode='<B>'+boxThresholdLine.name+':</B><INPUT type="text" id="'+thresholdInputID+'" size="5" value="'+boxThresholdLine.initValue+'"><INPUT type="button" onclick="idvis.registeredCharts['+BP.chartID+'].changeThreshold(document.getElementById(\''+thresholdInputID+'\').value)" value="Apply">';
			if (bpData.threshold.action) {
				htlmCode+='&nbsp;&nbsp;<INPUT type="button" value="'+bpData.threshold.action[0]+'" onclick="idvis.registeredCharts['+BP.chartID+'].getExcludedBoxes()">';
			}
			document.getElementById(formDivID).innerHTML=htlmCode;
		}

		/* Canvas chart*/
		this.canvas=Raphael(canvasDivID,canvasWidth,canvasHeight);

		/* Back panel */
		this.canvas.rect(0,0,canvasWidth,canvasHeight,0).attr({fill:'#fff',stroke:'#000'});

		/* Chart area */
		this.canvas.rect(boxPlotX0,boxPlotY0,boxPlotW,boxPlotH,0).attr({fill:'#F3F3F3',stroke:'#000'}); //,cursor:"crosshair"

		/* highlight zone */
		highlightZone=(horizontalBP)? this.canvas.rect(boxPlotX0,boxPlotY0,boxPlotW,boxSpace,0) : this.canvas.rect(boxPlotX0,boxPlotY0,boxSpace,boxPlotH,0);
		highlightZone.attr({fill:'#F00',stroke:'none','fill-opacity':0.3}).hide();

		/* Display chart */
		let boxPlotSize=(horizontalBP)? boxPlotW : boxPlotH,
			failed=false,failedValue=null;
		if (!idvis.isNumber(toScale(minBoxPlotValue))) {failed=true; failedValue=minBoxPlotValue;}
		else if (!idvis.isNumber(toScale(maxBoxPlotValue))) {failed=true; failedValue=maxBoxPlotValue;}
		if (failed) {
			alert('ERROR: Conversion of '+failedValue+' is not a valid number.');
			this.canvas.text(boxPlotX0+(boxPlotW/2),boxPlotY0+(boxPlotH/2),'Drawing has failed.').attr({'font-weight':'bold','font-size':14,fill:'#DD0000'});
			return;
		}
		let optimize=(this.forceTo0==1)? false : true,
			usedMinValue=(this.forceTo0)? 0 : toScale(minBoxPlotValue),
			optRange=idvis.lib.getChartScaleRange(optimize,usedMinValue,toScale(maxBoxPlotValue),toScale(boxPlotSize)), // requires chartLibrary
			optMinV=optRange[1],
			tickSizeV=optRange[3],
			boxPlotStartY0=boxPlotY0+boxPlotH; // verticalBP only
		minValue=optRange[0];
		maxValue=optRange[2];
		pix2valScale=(maxValue-minValue)/boxPlotSize;

		/** Axes titles & ticks **/
		/* Titles */
		let axisXtext,axisYtext;
		if (horizontalBP) {
			axisXtext=(bpData.valuesTitle)? bpData.valuesTitle : 'Values';
			axisYtext=(bpData.boxesTitle)? bpData.boxesTitle : 'Boxes';
		}
		else {
			axisYtext=(bpData.valuesTitle)? bpData.valuesTitle : 'Values';
			axisXtext=(bpData.boxesTitle)? bpData.boxesTitle : 'Boxes';
		}
		this.canvas.text(boxPlotX0+boxPlotW/2,canvasHeight-15,axisXtext).attr({'font-weight':'bold','font-size':14});
		let ty=boxPlotY0+(boxPlotH/2);
		this.canvas.text(15,ty,axisYtext).attr({'font-weight':'bold','font-size':14}).rotate(-90,15,ty);
		/* ticks */
		let posX,posY,
			fixed=(tickSizeV < 1)? (Math.round(1/tickSizeV)+'').length : (Math.round(tickSizeV)==tickSizeV)? 0 : 1;
		if (horizontalBP) {// X ticks
			posY=boxPlotStartY0+1;
			posX=boxPlotX0;
			let tickV=optMinV;
			while (tickV <= maxValue) {
				posX=boxPlotX0+Math.round((tickV-minValue)/pix2valScale);
				if (posX >= boxPlotX0) {
					this.canvas.path('M'+posX+' '+posY+' l0 5');
					this.canvas.text(posX,posY+10,tickV.toFixed(fixed));
				}
				tickV+=tickSizeV;
			}
		}
		else { // Y axis
			posX=boxPlotX0-6;
			let tickV=optMinV;
			while (tickV <= maxValue) {
				posY=boxPlotStartY0-Math.round((tickV-minValue)/pix2valScale);
				if (posY <= boxPlotStartY0) {
					this.canvas.path('M'+posX+' '+posY+' l5 0');
					this.canvas.text(posX-1,posY,tickV.toFixed(fixed)).attr({'text-anchor':'end'});
				}
				tickV+=tickSizeV;
			}
		}
//console.log('Canvas='+(new Date().getTime()-t0));

		/***** Drawing boxes & labels *****/
		let outlierRadius=Math.min(3,Math.round(boxSize/2)),
			whiskersW=Math.max(7,Math.round(boxSize/2)),
			boxPos=(horizontalBP)? boxPlotY0-(boxSpace/2) : boxPlotX0-(boxSpace/2);
		for (let i=0; i<boxList.length; i++) {
			boxPos+=boxSpace;
			/* label */
			let lx,ly,la;
			if (horizontalBP) {
				lx=boxPlotX0-5;
				ly=boxPos;
				la=0;
			}
			else {
				lx=boxPos;
				ly=boxPlotStartY0+5;
				la=-90; //-45;
			}
			boxList[i].labelSvg=BP.canvas.text(lx,ly,boxList[i].label)
			.attr({'font-size':labelFontSize,'font-weight':'bold','text-anchor':'end','fill':boxList[i].color})
			.rotate(la,lx,ly).data('js',boxList[i])
			.hover(function(){setBoxEmphasis(this.data('js'),'on');},function(){setBoxEmphasis(this.data('js'),'off');})
			.click(function() {boxClicked(this.data('js'));});
			//if (boxOnClick) {boxList[i].labelSvg.click(function() {boxOnClick(this.data('js').id)});}

			/* box */
			let box,median,whiskers,
				lowerWhSize=Math.round((toScale(boxList[i].lowerQuartile)-toScale(boxList[i].lowerWhisker))/pix2valScale),
				upperWhSize=Math.round((toScale(boxList[i].upperWhisker)-toScale(boxList[i].upperQuartile))/pix2valScale);
			if (horizontalBP) {
				let bx=boxPlotX0+Math.round((toScale(boxList[i].lowerQuartile)-minValue)/pix2valScale),
					by=Math.round(boxPos-(boxSize/2)),
					bw=Math.round((toScale(boxList[i].upperQuartile)-toScale(boxList[i].lowerQuartile))/pix2valScale);
				box=BP.canvas.rect(bx,by,bw,boxSize);
				let medx=boxPlotX0+Math.round((toScale(boxList[i].median)-minValue)/pix2valScale);
				median=BP.canvas.path('M'+medx+' '+by+' l0 '+boxSize);
				let wlx=boxPlotX0+Math.round((toScale(boxList[i].lowerWhisker)-minValue)/pix2valScale),
					wy=Math.round(boxPos-(whiskersW/2)),
					wux=boxPlotX0+Math.round((toScale(boxList[i].upperWhisker)-minValue)/pix2valScale);
				whiskers=BP.canvas.path('M'+wlx+' '+wy+' l0 '+whiskersW+' M'+wlx+' '+boxPos+' l'+lowerWhSize+' 0 M'+wux+' '+wy+' l 0 '+whiskersW+' M'+wux+' '+boxPos+' l-'+upperWhSize+' 0');
				boxList[i].x=wux;
				boxList[i].y=boxPos;
				let meax=boxPlotX0+Math.round((toScale(boxList[i].mean)-minValue)/pix2valScale);
				boxList[i].meanSvg=BP.canvas.path('M'+meax+' '+by+' l0 '+boxSize);
			}
			else {
				let bx=Math.round(boxPos-(boxSize/2)),
					by=boxPlotStartY0-Math.round((toScale(boxList[i].upperQuartile)-minValue)/pix2valScale),
					bh=Math.round((toScale(boxList[i].upperQuartile)-toScale(boxList[i].lowerQuartile))/pix2valScale);
				box=this.canvas.rect(bx,by,boxSize,bh);
				let medy=boxPlotStartY0-Math.round((toScale(boxList[i].median)-minValue)/pix2valScale);
				median=BP.canvas.path('M'+bx+' '+medy+' l'+boxSize+' 0');
				let wx=Math.round(boxPos-(whiskersW/2)),
					wly=boxPlotStartY0-Math.round((toScale(boxList[i].lowerWhisker)-minValue)/pix2valScale),
					wuy=boxPlotStartY0-Math.round((toScale(boxList[i].upperWhisker)-minValue)/pix2valScale);
				whiskers=BP.canvas.path('M'+wx+' '+wly+' l'+whiskersW+' 0 M'+boxPos+' '+wly+' l0 -'+lowerWhSize+' M'+wx+' '+wuy+' l'+whiskersW+' 0 M'+boxPos+' '+wuy+' l0 '+upperWhSize);
				boxList[i].x=boxPos;
				boxList[i].y=wuy;
				let meay=boxPlotStartY0-Math.round((toScale(boxList[i].mean)-minValue)/pix2valScale);
				boxList[i].meanSvg=BP.canvas.path('M'+bx+' '+meay+' l'+boxSize+' 0');
			}
			//box.attr({fill:boxList[i].color,stroke:boxList[i].color,'fill-opacity':0.3});
			median.attr({'stroke-width':3}); //stroke:boxList[i].color,
			boxList[i].meanSvg.attr({'stroke-width':3,'stroke-dasharray':'.',stroke:boxList[i].color});
			//whiskers.attr({stroke:boxList[i].color,'stroke-width':2});
			boxList[i].svg=BP.canvas.set(box,median,whiskers).attr({fill:boxList[i].color,stroke:boxList[i].color,'fill-opacity':0.5})
			.data('js',boxList[i])
			.hover(function(){setBoxEmphasis(this.data('js'),'on');},function(){setBoxEmphasis(this.data('js'),'off');}) // this=any element in set (not set!)
			.click(function() {boxClicked(this.data('js'));});
			if (boxList[i].excluded) {
				boxList[i].svg.attr({'stroke-dasharray':'- '}); //,'stroke-opacity':0.5
				boxList[i].meanSvg.attr({'stroke-opacity':0.5});
			}
			//if (boxOnClick) {boxList[i].svg.click(function() {boxOnClick(this.data('js').id)});}

			/* outliers */
			for (let j=0; j< boxList[i].outlierList.length; j++) {
				let outlier=boxList[i].outlierList[j];
				if (horizontalBP) {
					let ox=boxPlotX0+Math.round((toScale(outlier.value)-minValue)/pix2valScale);
					outlier.svg=BP.canvas.circle(ox,boxPos,outlierRadius);
					outlier.x=ox;
					outlier.y=boxPos;
				}
				else {
					let oy=boxPlotStartY0-Math.round((toScale(outlier.value)-minValue)/pix2valScale);
					outlier.svg=BP.canvas.circle(boxPos,oy,outlierRadius);
					outlier.x=boxPos;
					outlier.y=oy;
				}
				outlier.svg.attr({fill:boxList[i].color,stroke:boxList[i].color,'fill-opacity':0.3})
				.data('js',outlier)
				.hover(function(){setOutlierEmphasis(this.data('js'),'on');},function(){setOutlierEmphasis(this.data('js'),'off');})
				.click(function(){selectOutlier(this.data('js'));}); //,'auto'
				if (boxList[i].excluded) {outlier.svg.attr({'stroke-dasharray':'. ','stroke-opacity':0.5});}
			}
		}

		/***** Threshold line *****/
		if (boxThresholdLine) {
			thresholdValue=boxThresholdLine.initValue;
			if (boxThresholdLine.axis=='X') { // horizontal BP
				boxThresholdLine.pathPos=boxPlotX0+Math.round((toScale(boxThresholdLine.initValue)-minValue)/pix2valScale);
				boxThresholdLine.path=BP.canvas.path('M'+boxThresholdLine.pathPos+' '+boxPlotY0+' l0 '+boxPlotH);
				if (boxThresholdLine.label && !editableBP) {
					boxThresholdLine.path.hover(
						function(){let th=this.data('js'); th.popup=idvis.lib.drawLabel(BP.canvas,this,th.pathPos,boxPlotStartY0,1,1,th.label,th.color);}, // chartLibrary
						function(){let th=this.data('js'); th.popup.remove(); th.popup=null;}
					);
				}
			}
			else { //Y axis (vertical BP)
				boxThresholdLine.pathPos=boxPlotStartY0-Math.round((toScale(boxThresholdLine.initValue)-minValue)/pix2valScale);
				boxThresholdLine.path=BP.canvas.path('M'+boxPlotX0+' '+boxThresholdLine.pathPos+' l'+boxPlotW+' 0');
				if (boxThresholdLine.label && !editableBP) {
					boxThresholdLine.path.hover(
						function(){var th=this.data('js'); th.popup=idvis.lib.drawLabel(BP.canvas,boxPlotX0,th.pathPos,1,th.label,th.color);}, // chartLibrary
						function(){var th=this.data('js'); th.popup.remove(); th.popup=null;}
					);
				}
			}
			boxThresholdLine.path.attr({stroke:boxThresholdLine.color,'stroke-width':2,'stroke-dasharray':boxThresholdLine.dash})
			.data({js:boxThresholdLine});
			boxThresholdLine.prevPos=boxThresholdLine.pathPos;
		}
	};

	/*********************** Events management ******************************/
	function setBoxEmphasis(boxObj,action) {
		let color=(action=='on')? '#F00' : boxObj.color;
		boxObj.labelSvg.attr({fill:color});
		//boxObj.svg.attr({fill:color,stroke:color});
		//for (var i=0; i< boxObj.outlierList.length; i++) {
		//	if (!boxObj.outlierList[i].isSelected) {boxObj.outlierList[i].svg.attr({fill:color,stroke:color});}
		//}
		if (action=='on') {
			if (horizontalBP) {highlightZone.attr({y:boxPlotY0+(boxSpace*boxObj.index)});} else {highlightZone.attr({x:boxPlotX0+(boxSpace*boxObj.index)});}
			highlightZone.show();
			if (boxObj.popupSvg) {boxObj.popupSvg.show();}
			else {
				let text=boxObj.label+':\n'+boxObj.numValues+' values\nUpper W='+boxObj.upperWhisker+'\nUpper Q='+boxObj.upperQuartile+'\nMedian='+boxObj.median+'\nLower Q='+boxObj.lowerQuartile+'\nLower W='+boxObj.lowerWhisker+'\nMean='+boxObj.mean.toPrecision(4)+'\n'+boxObj.outlierList.length+' outlier(s)';
				boxObj.popupSvg=idvis.lib.drawLabel(BP.canvas,boxObj.x,boxObj.y,1,text,boxObj.color);
			}
		}
		else {
			boxObj.popupSvg.hide();
			highlightZone.hide();
		}
	}

	function boxClicked(boxObj) {
		if (idvis.modKeyPressed && excludeOnClick) {
			if (excludeOnClick(boxObj.id,boxObj.excluded)) { // 2 parameters boxID,isExcluded
				let update=(boxObj.excluded)? [false,'',''] : [true,'- ','. '];
				boxObj.excluded=update[0];
				boxObj.svg.attr({'stroke-dasharray':update[1]});
				for (let j=0; j<boxObj.outlierList.length; j++) {
					boxObj.outlierList[j].svg.attr({'stroke-dasharray':update[2]});
				}
			}
			else {alert('ERROR: External resource indicates update failure!');}
		}
		else if (boxOnClick) {boxOnClick(boxObj.id);}
	}

	this.changeThreshold=function(newValue) {
		let scaledNewValue=toScale(newValue);
		if (!idvis.isNumber(scaledNewValue)) { // chartLibrary
			alert("'"+newValue+"' is not a valid number!");
			document.getElementById(thresholdInputID).value=thresholdValue;
			return;
		}
		if (scaledNewValue < minValue || scaledNewValue > maxValue) {
			alert("'"+newValue+"' is out of range!");
			document.getElementById(thresholdInputID).value=thresholdValue;
			return;
		}
		boxThresholdLine.setValues(newValue);
		thresholdValue=newValue;
		//var shift=Math.round((toScale(newValue)-toScale(boxThresholdLine.startValue))/pix2valScale);
		if (boxThresholdLine.axis=='X') { // horizontal BP
			let shift=Math.round((scaledNewValue-toScale(boxThresholdLine.startValue))/pix2valScale);
			boxThresholdLine.path.transform('t'+shift+',0');
		}
		else { // vertical BP
			let shift=Math.round((toScale(boxThresholdLine.startValue)-scaledNewValue)/pix2valScale);
			boxThresholdLine.path.transform('t0,'+shift);
		}
	};

	this.getExcludedBoxes=function() {
		let excludedBoxes=[];
		for (let i=0; i<boxList.length; i++) {
			if (boxList[i].median < boxThresholdLine.initValue) excludedBoxes.push(boxList[i].id);
		}
		if (!bpData.threshold.action[1](boxThresholdLine.initValue,excludedBoxes)) { // call to external function must return a non-0/null value
			alert('ERROR: External resource indicates update failure!');
			return;
		}
		// OK to exclude
		for (let i=0; i<boxList.length; i++) {
			let update=[];
			if (boxList[i].median < boxThresholdLine.initValue) {
				if (!boxList[i].excluded) {
					update=[true,'- ',0.5,'. '];
				}
			}
			else if (boxList[i].excluded) { // median >= threshold
				update=[false,'',1,''];
			}
			if (update.length) {
				boxList[i].excluded=update[0];
				boxList[i].svg.attr({'stroke-dasharray':update[1],'stroke-opacity':update[2]});
				for (let j=0; j<boxList[i].outlierList.length; j++) {
					boxList[i].outlierList[j].svg.attr({'stroke-dasharray':update[3],'stroke-opacity':update[2]});
				}
			}
		}
	};

	function setOutlierEmphasis(outObj,action) {
		if (action=='on') {
			if (outObj.popupSvg) {outObj.popupSvg.show();}
			else {
				let text=(outObj.label)? outObj.label+': ' : '';
				text+=outObj.value;
				outObj.popupSvg=idvis.lib.drawLabel(BP.canvas,outObj.x,outObj.y,4,text,outObj.box.color); // chartLibrary
			}
			outObj.svg.attr({fill:'#F00',stroke:'#F00'});
		}
		else if (!outObj.isSelected) { // off
			outObj.popupSvg.hide();
			outObj.svg.attr({fill:outObj.box.color,stroke:outObj.box.color});
		}
	}

	function selectOutlier(outObj) { //,action
		//if (action=='auto') { // invert selection
			if (outObj.isSelected) { // already showing label => set OFF
				outObj.isSelected=false;
				setOutlierEmphasis(outObj,'off');
			}
			else { // set ON
				outObj.isSelected=true;
				setOutlierEmphasis(outObj,'on');
			}
		//}
/*
		else if (action=='on') { // set ON
			if (!outObj.isSelected) {
				outObj.isSelected=true;
				setOutlierEmphasis(outObj,'on');
			}
		}
		else { // OFF
			outObj.isSelected=false;
			setOutlierEmphasis(outObj,'off');
		}
*/
	}

/*================================ Nested objects =============================*/

	/*****************************************************************
					 Outlier point object
	******************************************************************/
	function outlierPoint(box,pointData,type) {
		this.box=box;
		this.value=pointData[0];
		this.label=(pointData[1])? pointData[1] : null;
		this.id=(pointData[2])? pointData[2] : this.label;
		this.type=type; // - for lower, + for upper
		this.x=this.y=null;
		this.svg=null; // graphical object
		this.popupSvg=null;
		this.isSelected=false;
	}

	/*****************************************************************
					  Box object
	******************************************************************/
	function box(BP,boxIdx,boxInfo) {
		this.index=boxIdx;
		this.label=boxInfo.label;
		this.id=(boxInfo.id)? boxInfo.id : boxIdx;
		this.excluded=(boxInfo.excluded)? true : false;
		this.color=(boxInfo.color)? boxInfo.color : '#000';
		this.numValues=null;
		this.median=null;
		this.mean=null;
		this.lowerQuartile=null;
		this.upperQuartile=null;
		this.lowerWhisker=null;
		this.upperWhisker=null;
		this.extremeValues=[];
		this.outlierList=[];
		this.computeBox(BP,boxInfo);
		this.x=this.y=null;
		this.labelSvg=null; // graphical label
		this.svg=null; // set of graphical objects
		this.meanSvg=null;
		this.popupSvg=null;
	//console.log(this);
	}
	box.prototype = {
		computeBox: function(BP,boxInfo) {
			let tmp2dArray=[],
				numErrors=0;
			for (let i=0; i<boxInfo.data.length; i++) { // converting into a 2D array
				let arr=(boxInfo.data[i]+'').split(':'); // to string (fails if true number)
				arr[0]*=1;
				if (idvis.isNumber(arr[0])) {tmp2dArray.push(arr);}
				else {numErrors++;}
			}
			if (numErrors>0) {alert('WARNING: '+numErrors+' non-valid value(s) found for '+this.label);}
			this.numValues=tmp2dArray.length;
			let dataList=tmp2dArray.sort(function(a,b){return a[0]-b[0];});
			this.extremeValues=[dataList[0][0],dataList[dataList.length-1][0]];
			//Box
			let medianRk=this.numValues/2, minMedianRk=Math.floor(medianRk), maxMedianRk=Math.ceil(medianRk);
			this.median=(medianRk != minMedianRk)? dataList[minMedianRk][0] : (dataList[medianRk-1][0]+dataList[medianRk][0])/2;
			let boxMinRk=maxMedianRk/2, minBoxMinRk=Math.floor(boxMinRk);
			this.lowerQuartile=(boxMinRk != minBoxMinRk)? dataList[minBoxMinRk][0] : (dataList[boxMinRk-1][0]+dataList[boxMinRk][0])/2;
			let boxMaxRk=this.numValues-boxMinRk, minBoxMaxRk=Math.floor(boxMaxRk);
			this.upperQuartile=(boxMaxRk != minBoxMaxRk)? dataList[minBoxMaxRk][0] : (dataList[boxMaxRk-1][0]+dataList[boxMaxRk][0])/2;
			this.mean=0;
			for (let i=0; i<this.numValues; i++) {this.mean+=dataList[i][0];}
			this.mean/=this.numValues;
//console.log(dataList.join(' '));
//console.log(medianRk+':'+minMedianRk+', '+boxMinRk+':'+minBoxMinRk+', '+boxMaxRk+':'+minBoxMaxRk);
			//Whiskers
			if (BP.outlierRule=='IQR') {
				let maxWhkSize=(this.upperQuartile-this.lowerQuartile)*1.5;
				//Lower
				let minWhkValue=this.lowerQuartile-maxWhkSize;
				for (let i=0; i<this.numValues; i++) {
					if (dataList[i][0] >= minWhkValue) {
						this.lowerWhisker=dataList[i][0];
						break;
					}
					else { // outlier
						this.outlierList.push(new outlierPoint(this,dataList[i],'-'));
					}
				}
				//Upper
				let maxWhkValue=this.upperQuartile+maxWhkSize;
				for (let i=this.numValues-1; i>=0; i--) {
					if (dataList[i][0] <= maxWhkValue) {
						this.upperWhisker=dataList[i][0];
						break;
					}
					else { // outlier
						this.outlierList.push(new outlierPoint(this,dataList[i],'+'));
					}
				}
			}
			else { // eg. 0.09 or 0.02
				//Lower
				let firstWhkIdx=Math.ceil(this.numValues*BP.outlierRule)-1;
//console.log(firstWhkIdx);
				this.lowerWhisker=dataList[firstWhkIdx][0];
				for (let i=0; i<firstWhkIdx; i++) {this.outlierList.push(new outlierPoint(this,dataList[i],'-'));}
				//upper
				let lastWhkIdx=this.numValues-firstWhkIdx-1;
//console.log(lastWhkIdx);
				this.upperWhisker=dataList[lastWhkIdx][0];
				for (let i=this.numValues-1; i>lastWhkIdx; i--) {this.outlierList.push(new outlierPoint(this,dataList[i],'-'));}
			}
		}
	};

};

/*
####>Revision history<####
# 2.0.3 Support for form positioning (PP 05/05/21)
# 2.0.2 'owner' svg data attribute replaced by 'js' (09/12/19)
# 2.0.1 [fix] Missing parameters in drawLabel call (PP 02/10/17)
# 2.0.0 Updated boxPlot.js to code-encapsulated idvis.boxPlot.js (PP 31/05/17)
# 1.1.0 Added forceTo0 option (PP 24/03/17)
# 1.0.9 Added this.getDivID function (PP 12/09/16)
# 1.0.8 Better handling of boxe size for user-defined chart width (PP PP 16/06/16)
# 1.0.7 Mean is also displayed (PP 13/04/16)
# 1.0.6 Uses chart registering for dynamic html calls (PP 20/02/16)
# 1.0.5 GPL license (PP 23/09/13)
# 1.0.4 Better label width calculation & excluded boxes display (PP 11/03/13)
# 1.0.3 Valid number checks on raw & toScale() numbers (PP 13/02/13)
# 1.0.2 Handles box exclusion/inclusion on shift+click (PP 05/02/13)
# 1.0.1 Added box visual exclusion & click event (PP 11/01/13)
# 1.0.0 First working version (PP 03/01/12)
*/