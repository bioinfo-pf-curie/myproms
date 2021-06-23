/*
################################################################################
# idvis.heatMap.js        3.2.0                                                #
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
                  HeatMap object
******************************************************************/
idvis.heatMap = function(mapData) {
    var HM=this;
    this.chartType='heatMap';
	this.chartID=idvis.lib.registerChart(this);

    /*=========== Configuration variables ==========*/
    var mainDivID=mapData.div,
        mainDiv=document.getElementById(mainDivID),
		paperDivID=mainDivID+'_paper',
		canvasID=mainDivID+'_canvas',
		formDivID;
	this.exportAsImage=mapData.exportAsImage;
	this.editMenu=null;
    if (mapData.editable || mapData.editableItems || mapData.searchable || this.exportAsImage) { // null if not declared
		this.editMenu=(mapData.editableItems)? mapData.editableItems : (mapData.editable)? {scope:true, type:false, color:true} : null; // used by external JS to know is form exists
		formDivID=mainDivID+'_form';
    }

	this.getDivID=function(){return paperDivID;};
    this.normProcess={ // default
		scope:'row', // 'row' 'column' 'all'
		reference:'mean', // z-score,median,mean,mid-range,min,max,user
		colors:'GbR' // GbR,GYR,GWR,BbR,BYR,BWR,BGY...
    };
	var maxRange=mapData.maxRange || {}, // max range for cell values {type:<absolute|quartile>,value:<[min,max]|number of quartiles>}
		rangeMinValue,rangeMaxValue;
	if (mapData.entities) { // Which entities are reprensented by row/column
		this.entities=mapData.entities;
		if (!this.entities.row) this.entities.row='Row';
		if (!this.entities.column) this.entities.row='Column';
	}
	else {this.entities={row:'Row',column:'Column'};}

	/* Heatmap colors management */
	const colorPatterns=[
						 ['GbR','green/black/red'],['GYR','green/yellow/red'],['GWR','green/white/red'],['GbB','green/black/blue'],
						 ['RbG','red/black/green'],['RbB','red/black/blue'],
						 ['BbR','blue/black/red'],['BYR','blue/yellow/red'],['BWR','blue/white/red'],['BbG','blue/black/green'],
						 ['BGY','blue/green/yellow']
						];
	const colorCharCode={R:'#F00',G:'#0F0',B:'#00F',b:'#000',Y:'#FF0',W:'#FFF'}; // B: blue, b: black
	const colorPairSet={ // not the same DOWN and UP
		DOWN:{
			Gb:'rgb(0%,+%,0%)',
			GY:'rgb(-%,100%,0%)', //!
			GW:'rgb(-%,100%,-%)',
			Bb:'rgb(0%,0%,+%)',
			BY:'rgb(-%,-%,+%)',
			BW:'rgb(-%,-%,100%)',
			BG:'rgb(0%,-%,+%)',
			Rb:'rgb(+%,0%,0%)'
		},
		UP:{
			bR:'rgb(+%,0%,0%)',
			YR:'rgb(100%,-%,0%)',
			WR:'rgb(100%,-%,-%)',
			GY:'rgb(+%,100%,0%)', //!
			bG:'rgb(0%,+%,0%)',
			bB:'rgb(0%,0%,+%)'
		}
	};

    if (mapData.normalization) {
		//for (let attr in mapData.normalization) {this.normProcess[attr] = mapData.normalization[attr];}
		Object.keys(mapData.normalization).forEach(function(attr) {HM.normProcess[attr] = mapData.normalization[attr];});
        if (!this.normProcess.scope.match('row|column|all')) {this.normProcess.scope='row';}
		if (!this.normProcess.reference.match('mean|median|mid-range|min|zero|user') || (this.normProcess.reference=='user' && !this.normProcess.limitValues)) {this.normProcess.reference='z-score';}
		//if (!this.normProcess.colors.match('GbR|GYR|GWR|BbR|BYR|BWR')) {this.normProcess.colors='GbR';} // B: blue, b: black
		let patternMatched=false;
		for (let i=0; i<colorPatterns.length; i++) {
			if (this.normProcess.colors==colorPatterns[i][0]) {
				patternMatched=true;
				break;
			}
		}
		if (!patternMatched) {this.normProcess.colors='GbR';}
		if (this.normProcess.reference=='z-score') {this.editMenu.scope=false;}
		//else if (this.normProcess.reference=='user') {this.editMenu.type=false;}
    }

	/* Cell value line */
	var cellsValueTrace=mapData.showTrace || false,
//this.treeSpace={row:0,column:0};
	/* Rows/columns layout & properties */
		rowLabelLocation=mapData.rowLabelPos || 'left',
		colLabelLocation=mapData.colLabelPos || 'top',
		colLabelOrientation=mapData.colLabelOrientation || 'vertical',
		colLabelAngle=(colLabelOrientation=='diagonal')? -45 : -90,
		colLabelRatio=Math.abs(Math.sin(colLabelAngle * (Math.PI / 180))), //radian
		movableLabel={},
		selectedLabels=[], // list of label IDs highlighted by clicking on a tree node
		minCellH=(mapData.noMinCellHeight)? 1 : 5;
		movableLabel.row=(mapData.moveLabel || mapData.moveRow)? true : false;
		movableLabel.column=(mapData.moveLabel || mapData.moveColumn)? true : false;

	/* Rows/columns font family */
	var rowLabelFontFamily,rowGroupLabelFontFamily,colLabelFontFamily,colGroupLabelFontFamily;
	if (mapData.fontFamily) {
		rowLabelFontFamily=mapData.fontFamily.rowLabels || mapData.fontFamily.labels || null;
		rowGroupLabelFontFamily=mapData.fontFamily.rowGroupLabels || mapData.fontFamily.labels || null;
		colLabelFontFamily=mapData.fontFamily.columnLabels || mapData.fontFamily.labels || null;
		colGroupLabelFontFamily=mapData.fontFamily.columnGroupLabels || mapData.fontFamily.labels || null;
	}

	/* Palette */
 	var palette={valueDistrib:[],trimValues:{left:0,right:0,svgLeft:null,svgRight:null},hasNegValues:false,hasPosValues:false,posX:null,posY:null,refLimit:null,svgElements:[]};
	//palette.posX=(rowLabelLocation=='left')?

	/* Row,column mouse events */
    var rowOnClick=mapData.rowOnClick || null,
		rowOnModClick=mapData.rowOnModClick || null,
		rowOnDblClick=mapData.rowOnDblClick || null,
		rowOnModDblClick=mapData.rowOnModDblClick || null,
		rowGroupOnClick=mapData.rowGroupOnClick || null,
		columnOnClick=mapData.columnOnClick || null,
		columnOnModClick=mapData.columnOnModClick || null,
		columnOnDblClick=mapData.columnOnDblClick || null,
		columnOnModDblClick=mapData.columnOnModDblClick || null,
		columnGroupOnClick=mapData.columnGroupOnClick || null,
		labelOnMove=mapData.labelOnMove || null; // if a row/column label was moved
	/* Cell mouse events */
	var oneCellSelection=(mapData.singleCellSelection)? true : false,
		lastSelectedCell,lastHighlightedCell, // js objects
		cellOnClick=mapData.cellOnClick || null,
		noCellExclusion=(mapData.noCellExclusion)? true : false,
		cellOnModClick=(!noCellExclusion && mapData.cellOnModClick)? mapData.cellOnModClick : null, // noCellExclusion takes over cellOnModClick
		cellOnMouseOver=mapData.cellOnMouseOver || null,
		annotationOnDelete=mapData.annotationOnDelete || null,
		//cellHighlightBgColor='#000', //(this.normProcess.colors.match('Y'))? '#FFF' : '#000';
		cellHighlightColor,noDataCellColor; // computed by computeCellsColor
	/* Tree mouse event */
	var labelsOnSelect=mapData.labelsOnSelect || {}; // {row/column:{text:...,action:...[,ranking:true/false]}}
	for (let type in {row:1,column:1}) {
		if (labelsOnSelect[type]) {
			if (!labelsOnSelect[type].action) {delete labelsOnSelect[type];}
			else if (!labelsOnSelect[type].text) {labelsOnSelect[type].text='Click for action';}
		}
	}
	/* Cell flagging */
	var flagText=mapData.flagText || 'Flagged',
		cellFlagVisible=true;
	var existFlaggedCells=false, // default. update in addRow()
		treePopup={}, // row/column popup SVG
		branchIsSelected=false, // Flag for tree branch selection
		branchColorIndex={row:{},column:{}},
		annotNameToColor={row:{},column:{}};
	//var annotColorIndex={row:{},column:{}};
	const colorList=['#000000','#0000FF','#4AA02C','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B']; // for clustering tree/label/annotation selection

    /*============= Data variables =============*/
    var columnList=[],rowList=[],groupList={row:[],column:[]};
	this.getColumnList=function() {return columnList;}; // to allow access from outside HM
	this.getRowList=function() {return rowList;}; // to allow access from outside HM
	var hasCells=false, // 1D or 2D
		branchList={}, // Clustering tree {row:[],column:[]}
		treeLength={},
		annotationList={row:{},column:{}},
		annotationCellSize={row:null,column:null},
		annotationLegends, // SVG set containing all SVG elements of legends
    //var numDefinedValues=0;
		maxColLabel='',maxRowLabel='',maxGroupAndLabels={column:'',row:''},
	//var maxColLabelSize=15; //7; // min.length used
	//var maxRowLabelSize=15; //7; // min.length used
		maxAbsZscore=0;

    /*============ Layout variables =============*/
	const topSpace=20,
		  bottomSpace=20,
		  leftSpace=20,
		  rightSpace=20;
	var paperWidth,paperHeight,hCellW,hCellH,flagSize,rowLabelFontSize,colLabelFontSize,rowLabelSpace,colLabelSpace,colDiagSpace,mapW,mapH,cellX0,cellY0,canvas,mapContext,
		midCell,defaultRadius,roundedCells, // for rounded cells
		paperWidth0,paperHeight0,background,backgroundPath,
		treeSpace={row:50,column:50},
		treeStart={row:0,column:0},
    //var mapGeometry={};
		labelPopupText,cellPopupText,
		hCellsSensor, //SVG rect over heatCells
		hCellFrame, // SVG frame around highlighted heatCell
		focusArea, //SVG drag rect over heatCells
		labelPosRanges={},
		maxOverlap={},
		dragParams={},
		dragContext=false,
		dblClickContext=false,
		helpIsVisible=false,
		helpLogo,helpPopup;

	/********** Set all columns ***************/
	this.setColumns=function(colData) {
		let rank=0;
		for (let i=0; i<colData.length; i++) {
			let newCol=new Label(this,'column',++rank,colData[i]);
			columnList.push(newCol);
			if (maxColLabel.length < newCol.label.length) maxColLabel=newCol.label;
		}
//console.log('COL: '+columnList.length);
	};

	/********** Add multiple or single row(s) & associated data ***************/
	this.addRows=function(rowsData) {
		for (let i=0; i<rowsData.length; i++) {
			HM.addRow(rowsData[i][0],rowsData[i][1],rowsData[i][2],rowsData[i][3]);
		}
	};
	this.addRow=function(rowInfo,rowData,excludedCells,flaggedCells) {
		let newRow=new Label(this,'row',rowList.length+1,rowInfo);
		rowList.push(newRow);
		if (maxRowLabel.length < newRow.label.length) maxRowLabel=newRow.label;
//console.log('ROW: '+rowData.length);
		for (let i=0; i<columnList.length; i++) { // columnList.length in case less values in current row
			let hCell=new HeatCell(this,rowData[i]);
			//hCell.row=newRow;
			hCell.rowIdx=rowList.length-1;
			//hCell.column=columnList[i];
			hCell.columnIdx=i;
			columnList[i].heatCellList.push(hCell);
			newRow.heatCellList.push(hCell);
			if (hCell.value !== null) {
				rangeMinValue=(rangeMinValue===undefined)? hCell.value : Math.min(rangeMinValue,hCell.value);
				rangeMaxValue=(rangeMaxValue===undefined)? hCell.value : Math.max(rangeMaxValue,hCell.value);
			}
		}
		// Excluded cells
		if (excludedCells) {
			for (let e=0; e<excludedCells.length; e++) {
				newRow.heatCellList[excludedCells[e]].isExcluded=true;
			}
		}
		// Flagged cells
		if (flaggedCells) {
			for (let e=0; e<flaggedCells.length; e++) {
				newRow.heatCellList[flaggedCells[e]].isFlagged=true;
				existFlaggedCells=true;
			}
		}
	};
	this.addRowsSecondValues=function(rowsData2) {
		for (let i=0; i<rowsData2.length; i++) {
			HM.addRowSecondValues(i,rowsData2[i]);
		}
	};
	this.addRowSecondValues=function(rowIdx,rowData2) {
		let row=rowList[rowIdx];
		for (let i=0; i<rowData2.length; i++) {
			row.heatCellList[i].value2=(idvis.isNumber(rowData2[i]))? rowData2[i] : null; // new attribute added to HeatCell
		}
	}

	/********** Define row or column groups of adjacent labels ***************/
	/* type: column or row
	 * groupData: array of objects[[label1,id1,popup,startIndex1,endIndex2],[...],...] */
	this.defineGroups=function(type,groupData) {
		movableLabel[type]=false; // disable label move
		//if (type=='column') {
		//	colLabelOrientation='vertical'; // too complicated for diagonal
		//	colLabelAngle=-90;
		//	colLabelRatio=1;
		//}
		if (!groupList[type]) groupList[type]=[];
		let labelList=(type=='row')? rowList : columnList;
		for (let i=0; i<groupData.length; i++) {
			let groupIdx=groupList[type].length,
				gr=new labelGroup(this,type,groupData[i]);
			groupList[type].push(gr);
			let maxLabelIdx=gr.labelIndexRange[0];
			for (let l=gr.labelIndexRange[0]; l<=gr.labelIndexRange[1]; l++) {
				labelList[l].groupIndex=groupIdx;
				if (labelList[maxLabelIdx].label.length < labelList[l].label.length) maxLabelIdx=l;
			}
			gr.labelIndexRange.push(maxLabelIdx);
			let maxLabel=gr.label+labelList[maxLabelIdx].label;
//console.log(maxLabelIdx,maxLabel);
			if (maxGroupAndLabels[type].length < maxLabel.length) maxGroupAndLabels[type]=maxLabel;
		}
	};


	/********** Add a row or column clustering tree ***************/
	this.addTree=function(type,dataStrg) {
		movableLabel[type]=false; // overwrites (just to be safe)
	    branchList[type]=[];
	    treeLength[type]=0;
	    let labelList=(type=='row')? rowList : columnList,
			labelIdList={}; // id -> Label
	    for (let i=0; i<labelList.length; i++) {
			labelIdList[labelList[i].id]=labelList[i];
	    }
	    let treeData=dataStrg.split(';');
	    for (let i=0; i<treeData.length; i++) {
			let branchData=treeData[i].split(','),
			//var leaf1=(branchData[1]>0)? branchList[type][branchData[1]-1] : labelIdList[-branchData[1]];
			//var leaf2=(branchData[2]>0)? branchList[type][branchData[2]-1] : labelIdList[-branchData[2]];
				leaf1=(branchData[1].match('-'))? labelIdList[branchData[1].substring(1)] : branchList[type][branchData[1]-1],
				leaf2=(branchData[2].match('-'))? labelIdList[branchData[2].substring(1)] : branchList[type][branchData[2]-1];
			//Checks if label ids are matched
			if (!leaf1 && branchData[1].match('-')) {alert('ERROR: id "'+branchData[1].substring(1)+'" not found in '+type+' labels while reading tree.'); break;}
			if (!leaf2 && branchData[2].match('-')) {alert('ERROR: id "'+branchData[2].substring(1)+'" not found in '+type+' labels while reading tree.'); break;}
			if (leaf1.rank > leaf2.rank) { // invert leaves
				let tmpLeaf=leaf1;
				leaf1=leaf2;
				leaf2=tmpLeaf;
			}
			//branchList[type].push(new TreeBranch(this,type,branchData[0],leaf1,leaf2,branchData[3]));
			branchList[type].push(new TreeBranch(this,type,leaf1.rank,leaf1,leaf2,branchData[3]));
			treeLength[type]=Math.max(treeLength[type],branchData[3]);
	    }
	};

	/********************* Heat map display (Raphael) *********************/
    this.draw=function() {
		if (branchList.row) {treeSpace.row=Math.min(500,Math.max(150,Math.round(5*branchList.row[branchList.row.length-1].distance/branchList.row[0].distance)));} /* [100 >= 5px for minDist <= 500] */
		if (branchList.column) {treeSpace.column=Math.min(500,Math.max(150,Math.round(5*branchList.column[branchList.column.length-1].distance/branchList.column[0].distance)));} /* [100 >= 5px for minDist <= 500] */
//console.log(treeSpace);
		let maxMapWidth=(treeSpace.row)? 1500-treeSpace.row : 1500;
		//Columns
		let numColumns=columnList.length;
		//hCellW=Math.round(maxMapWidth/numColumns);
		//if (hCellW < 11) {hCellW=11;} else if (hCellW > 50) {hCellW=50;}
		hCellW=(mapData.cellWidth)? mapData.cellWidth : (mapData.width)? mapData.width/numColumns : Math.max(11,Math.min(50,Math.round(maxMapWidth/numColumns))); // mapData.cellWidth overwrites mapData.width!!!
		//Rows
		let numRows=rowList.length;
		if (numRows) { // 2D
			hasCells=true;
			//hCellH=Math.round(750/numRows);
			//if (hCellH < minCellH) {hCellH=minCellH;} else if (hCellH > 25) {hCellH=25;} // 13
			//if (hCellH > hCellW) {hCellH=hCellW;}
			hCellH=(mapData.cellHeight)? mapData.cellHeight : (mapData.height)? mapData.height/numRows : Math.min(hCellW,Math.max(minCellH,Math.min(25,Math.round(750/numRows))));
			if (mapData.roundedCells) {
				hCellW=hCellH=Math.min(hCellW,hCellH); // same H & W
				roundedCells={range:[],fitMainValue:mapData.roundedCells.fitMainValue || false,sizeLabel:'Size'};
				if (mapData.roundedCells.range) {roundedCells.range=mapData.roundedCells.range;}
				if (mapData.roundedCells.sizeLabel) {roundedCells.sizeLabel=mapData.roundedCells.sizeLabel;}
				midCell=hCellW/2;
				defaultRadius=Math.round(midCell)-1.5;
			}
			maxOverlap.row=hCellH/2;
		}
		else {hCellH=0;} // 1D
		maxOverlap.column=hCellW/2;

		flagSize=(hCellH >= 5)? Math.min(10,Math.min(hCellW,hCellH)) : 3; //Size of flag triangle
//console.log('FLAG',flagSize);
//var rowTreeSpace=(branchList.row)? treeSpace : 0;
//var columnTreeSpace=(branchList.column)? treeSpace : 0;

		//Column labels space
		colLabelFontSize=(hCellW >= 25)? 14 : (hCellW >= 20)? 12 : 8;
		//Row labels space
		rowLabelFontSize=(hCellH >= 25)? 14 : (hCellH >= 20 || !hasCells)? 12 : 8; // 12 default size if 1D

		/* Temporary paper to compute label space */
		let tmpPaper=Raphael(mainDivID,0,0),
			colText=(maxGroupAndLabels.column && maxGroupAndLabels.column.length > maxColLabel.length)? maxGroupAndLabels.column : maxColLabel,
			colLabel=tmpPaper.text(0,0,colText).attr({'font-size':colLabelFontSize,'font-weight':'bold'});
		if (colLabelFontFamily) {colLabel.attr({'font-family':colLabelFontFamily});}
		colLabelSpace=colLabel.getBBox().width;
		if (colGroupLabelFontFamily) {
			colLabel.attr({'font-family':colGroupLabelFontFamily});
			colLabelSpace=Math.max(colLabelSpace,colLabel.getBBox().width);
		}
		if (maxGroupAndLabels.column) colLabelSpace+=15;
		if (colLabelOrientation=='vertical') {
			colDiagSpace=0;
		}
		else { // diagonal
			colLabelSpace=Math.round(colLabelSpace*colLabelRatio);
			//colDiagSpace=(!branchList.row && ((colLabelLocation=='top' && rowLabelLocation=='left') || (colLabelLocation=='bottom' && rowLabelLocation=='right')))? Math.round((colLabelSpace*Math.abs(Math.cos(colLabelAngle * (Math.PI / 180))))-(hCellW/2)) : 0;
			//colDiagSpace=((colLabelLocation=='top' && rowLabelLocation=='right') || (colLabelLocation=='bottom' && rowLabelLocation=='right'))? Math.round(colLabelSpace-(hCellW/2)) : 0;
			colDiagSpace=(rowLabelLocation=='right')? Math.round(colLabelSpace-(hCellW/2)) : 0;
		}
		colLabelSpace=Math.max(colLabelSpace,100); // room for palette
		if (hasCells) {
			let rowText=(maxGroupAndLabels.row && maxGroupAndLabels.row.length > maxRowLabel.length)? maxGroupAndLabels.row : maxRowLabel,
				rowLabel=tmpPaper.text(0,0,rowText).attr({'font-size':rowLabelFontSize,'font-weight':'bold'});
			if (rowLabelFontFamily) {rowLabel.attr({'font-family':rowLabelFontFamily});}
			rowLabelSpace=Math.max(rowLabel.getBBox().width,200); // room for palette
			if (rowGroupLabelFontFamily) {
				rowLabel.attr({'font-family':rowGroupLabelFontFamily});
				rowLabelSpace=Math.max(rowLabelSpace,rowLabel.getBBox().width);
			}
			if (maxGroupAndLabels.row) rowLabelSpace+=15;
		}
		else {rowLabelSpace=60;} // room for annotation labels
		tmpPaper.remove();

//console.log(hCellW+','+hCellH);

		annotationCellSize.row=Math.min(hCellW,20);
		annotationCellSize.column=(hasCells)? Math.max(10,Math.min(hCellH,20)) : 20; // 1/2D for type=column (same in deleteAnnotation)

		mapW=columnList.length*hCellW;
		mapH=rowList.length*hCellH;
		paperWidth=leftSpace+rowLabelSpace+mapW+treeSpace.row+colDiagSpace+rightSpace;
		paperHeight=topSpace+colLabelSpace+mapH+treeSpace.column+bottomSpace;
		paperWidth0=paperWidth; paperHeight0=paperHeight;
//console.log(paperWidth+','+paperHeight);
		cellX0=(rowLabelLocation=='left')? leftSpace+rowLabelSpace : (colLabelLocation=='bottom' && colLabelOrientation=='diagonal')? leftSpace+colDiagSpace+treeSpace.row : leftSpace+treeSpace.row;
		cellY0=(colLabelLocation=='top')? topSpace+colLabelSpace : topSpace+treeSpace.column;
//console.log(cellX0+','+cellY0);
		palette.posX=(rowLabelLocation=='left')? Math.max(cellX0-250,10) : Math.min(cellX0+mapW+colDiagSpace+40,paperWidth-220);
		palette.posY=(colLabelLocation=='top')? Math.max(cellY0-150,5) : Math.min(cellY0+mapH+50,paperHeight-120);

		/* DIVs */
		let mainDiv=document.getElementById(mainDivID),
			formDivHeight=0;
		if ((this.editMenu && hasCells) || this.exportAsImage) { // hasCells to check for 2D
			mainDiv.innerHTML='<DIV id="'+formDivID+'" style="margin-bottom:4px"></DIV>';
			/* Form menu */
			initializeHeatMapForm();
			formDivHeight=document.getElementById(formDivID).offsetHeight;
		}
		let divStrg='<DIV style="position:relative;">';
		divStrg+='<CANVAS id="'+canvasID+'" style="position:absolute;top:'+cellY0+'px;left:'+cellX0+'px" width="'+mapW+'" height="'+mapH+'"></CANVAS><DIV id="'+paperDivID+'" style="background:rgba(255,255,255,0);"></DIV>'; // class="transparent"
//divStrg+='<CANVAS id="'+canvasID+'_zoom" style="position:absolute;background-color:#0000FF;top:'+cellY0+'px;left:'+cellX0+'px" width="'+(mapW/2)+'" height="'+mapH+'"></CANVAS><DIV id="'+paperDivID+'_zoom" style="background:rgba(255,255,255,0);background-color:#00FFFF">Hello Patrick</DIV>'; // class="transparent"
		divStrg+='</DIV>';
		mainDiv.innerHTML+=divStrg;

		let paperDiv=document.getElementById(paperDivID);
		paperDiv.style.position='absolute'; // inside parent DIV (because relative)
		paperDiv.style.top='0px'; paperDiv.style.left='0px';

		/* svg & Canvas charts*/
		this.paper=Raphael(paperDivID,paperWidth,paperHeight);
		canvas = document.getElementById(canvasID);
		mapContext = canvas.getContext('2d');

		/* Adjust mainDiv size (because of position:relative of paperDiv) */
//console.log(formDivHeight, paperDiv.offsetWidth, paperDiv.offsetHeight);
		mainDiv.style.height=(formDivHeight + paperDiv.offsetHeight)+'px';
		mainDiv.style.width=paperDiv.offsetWidth+'px';

		/* Back panel */
		//this.paper.rect(0,0,paperWidth,paperHeight,0).attr({fill:'#fff',stroke:'#000'});
		//Relative
		//var backgroundPathStrg='M0,0 l'+paperWidth+',0 l0,'+paperHeight+',0 l-'+paperWidth+',0 l0,'+(cellY0-paperHeight)+' l'+cellX0+',0 l0,'+mapH+' l'+mapW+',0 l0,-'+mapH+' l-'+(mapW+cellX0)+',0 Z';
		//Absolute
		let backgroundPathStrg='M0,0 L'+paperWidth+',0 L'+paperWidth+','+paperHeight+' L0,'+paperHeight+' L0,'+cellY0+' L'+cellX0+','+cellY0+' L'+cellX0+','+(cellY0+mapH)+' L'+(cellX0+mapW)+','+(cellY0+mapH)+' L'+(cellX0+mapW)+','+cellY0+' L0,'+cellY0+' Z';
//console.log(backgroundPathStrg);
		backgroundPath=this.paper.path(backgroundPathStrg).attr({fill:'#fff',stroke:'none'});
//console.log(backgroundPath);
		background=this.paper.rect(0,0,paperWidth,paperHeight,0).attr({stroke:'#000'});

		/* Labels */
		drawLabels();

		/** 2D **/
		if (hasCells) {

			/* Correction for values max range */
			if (maxRange.type) {applyMaxRange();}

			/* Heatmap */
			this.computeHeatMap(); // also draw cells with color

			/* heatCell frame */
			hCellFrame=HM.paper.rect(1,1,hCellW-2,hCellH-2).attr({'stroke-width':2,stroke:cellHighlightColor}).hide(); // stroke color reset in computeCellsColor()

			/* Sensor */
			hCellsSensor=HM.paper.rect(cellX0,cellY0,mapW,mapH,0).attr({fill:'#555','fill-opacity':0,stroke:'none'})
			.mousemove(function(e){
							let hCell=findHeatCell(e);
							//if (!hCell.aspect) return; // mouseout too fast **********************************
							if (lastHighlightedCell) {
								//if (hCell.rowIdx==lastHighlightedCell.rowIdx && hCell.columnIdx==lastHighlightedCell.columnIdx) {return;} //same cell
								if (hCell===lastHighlightedCell) {return;} //same cell
								else {setCellEmphasis(lastHighlightedCell,'off');} // new cell
							}
				if (!roundedCells || hCell.value !== null) { // No cell highlighting for empty rounded cell
							setCellEmphasis(hCell,'on',e);
							lastHighlightedCell=hCell;
				}
						})
			.mouseout(function(){if (lastHighlightedCell) {setCellEmphasis(lastHighlightedCell,'off'); }})  // lastHighlightedCell=null; ( DO NOT SET TO null to allow pseudo-click event )
			.drag(function(dx,dy,x,y,e){dragFocusMove(dx,dy,x,y,e);}, //_extendDragging,
				  function(x,y,e){dragFocusStart(x,y,e);}, //_setDragging,
				  function(){dragFocusEnd();} //_endDragging
				  );
            /* !!! hCellsSensor.click handled at drag end level due to conflict in non-Firefox browsers !!!
			/*if (cellOnClick) {hCellsSensor.click(function(e) {if ((e.ctrlKey || e.shiftKey) && !noCellExclusion) {excludeCell(lastHighlightedCell);} else if (lastHighlightedCell) {selectCell(lastHighlightedCell);}});}
			/*else if (!noCellExclusion) {hCellsSensor.click(function(e){if (e.ctrlKey || e.shiftKey) {excludeCell(lastHighlightedCell);}});}
			*/

			/* Focus area & control(s) */
			let setY=7,
				closeSet=HM.paper.set( // Close icon
                    HM.paper.rect(8,setY,18,18,2).attr({stroke:'none'}),
                    HM.paper.path('M11,'+(setY+3)+' l12,12').attr({stroke:'#000','stroke-width':2}),
                    HM.paper.path('M23,'+(setY+3)+' l-12,12').attr({stroke:'#000','stroke-width':2})
                ).hide()
                .click(function() {closeFocusArea();});
			focusArea=HM.paper.rect(0,0,0,0,0)
			.attr({'stroke-width':2,'fill-opacity':0}) // 0.2
			.data({startX:0,startY:0,status:'off',labelIdx:[],
				  bgRow:HM.paper.rect(0,0,0,0,5).attr({fill:'#FFF',stroke:'none','fill-opacity':0.7}).hide(),
				  bgRow2:HM.paper.rect(0,0,0,0,5).attr({fill:'#FFF',stroke:'none','fill-opacity':0.7}).hide(),
				  bgColumn:HM.paper.rect(0,0,0,0,5).attr({fill:'#FFF',stroke:'none','fill-opacity':0.7}).hide(),
				  bgColumn2:HM.paper.rect(0,0,0,0,5).attr({fill:'#FFF',stroke:'none','fill-opacity':0.7}).hide(),
				  close:closeSet})
			//focusArea.dblclick(function() {zoomHeatMap()})
			.hide();
			//var setY=27;
			setY+=20;
			if (labelsOnSelect.row) {
				let selRowsSet=HM.paper.set( // Row list icon
					HM.paper.rect(8,setY,18,18,2).attr({stroke:'none'}), //.attr({fill:'#FF0'}),
					HM.paper.path('M11,'+(setY+4)+' l12,0').attr({stroke:'#000','stroke-width':2}),
					HM.paper.path('M11,'+(setY+9)+' l12,0').attr({stroke:'#000','stroke-width':2}),
					HM.paper.path('M11,'+(setY+14)+' l12,0').attr({stroke:'#000','stroke-width':2})
				).hide().
				click(function() {focussedLabelsAction('focus','row');});
				focusArea.data({selRows:selRowsSet});
				setY+=20;
			}
			if (labelsOnSelect.column) {
				let selColumnsSet=HM.paper.set( // Column list icon
					HM.paper.rect(8,setY,18,18,2).attr({stroke:'none'}), //.attr({fill:'#FF0'}),
					HM.paper.path('M12,'+(setY+3)+' l0,12').attr({stroke:'#000','stroke-width':2}),
					HM.paper.path('M17,'+(setY+3)+' l0,12').attr({stroke:'#000','stroke-width':2}),
					HM.paper.path('M22,'+(setY+3)+' l0,12').attr({stroke:'#000','stroke-width':2})
				).hide()
				.click(function() {focussedLabelsAction('focus','column');});
				focusArea.data({selColumns:selColumnsSet});
				setY+=20;
			}
		}

		/* Trees */
		treeStart.row=(rowLabelLocation=='left')? cellX0+mapW+2 : cellX0-2;
		treeStart.column=(colLabelLocation=='top')? cellY0+mapH+2 : cellY0-2;
		if (branchList.row) {drawTree('row');} // row tree
		if (branchList.column) {drawTree('column');}

		/* keyboard listener */
		//addKeyboardEvents(); handled by chartLibrary!

		/* Add color to cells */
		//this.computeHeatMap();

		/* Annotation legends */
		annotationLegends=HM.paper.set();

		/* Help popup */
		let helpY=(colLabelLocation=='top')? paperHeight-12 : 12,
			helpLink=this.paper.text(6,helpY,'HELP').attr({'text-anchor':'start','font-size':12,'font-weight':'bold',fill:'#FFF'}),
		//helpLink.click(showHideHelpText);
			hlBox=helpLink.getBBox(),
			hlb=this.paper.rect(hlBox.x-2,hlBox.y-1,hlBox.width+4,hlBox.height+2,0).attr({stroke:'none',fill:'#444','fill-opacity':1});
		//hlb.click(showHideHelpText);
		helpLogo=this.paper.set(helpLink,hlb);
		helpLink.toFront();
		this.paper.set(helpLink,hlb)
		.hover(function(){this.attr({cursor:'pointer'});},function(){this.attr({cursor:'auto'});})
		.click(showHideHelpText);
		let yShift=50;
		//var helpText='HELP:\n[Mouseover+/-Shift] for cell information.';
		var helpText='HELP:';
		if (hasCells) {
			helpText+='\n[Mouseover] for cell information.';
			if (!noCellExclusion) {helpText+='\n[Shift+Click] to exclude/include a cell.';}
			if (cellOnClick) {
				helpText+='\n[Click] to select/unselect a cell.';
				yShift+=7;
			}
			if (rowOnClick || rowOnModClick || rowOnDblClick || rowOnModDblClick) {
				helpText+='\n[Shift+/-Double Click] on '+HM.entities.row+' label for special action.';
				yShift+=7;
			}
		}
		if (columnOnClick || columnOnModClick || columnOnDblClick || columnOnModDblClick) {
			helpText+='\n[Shift+/-Double Click] on '+HM.entities.column+' label for special action.';
			yShift+=7;
		}
		if (movableLabel.row || movableLabel.column) {
			helpText+='\nGrab label to move ';
			helpText+=(movableLabel.row && movableLabel.column)? HM.entities.row+'/'+HM.entities.column : (movableLabel.row)? HM.entities.row : HM.entities.column;
			helpText+='.';
			yShift+=7;
		}
		if (branchList.row || branchList.column) {
			helpText+='\n[Click] on tree nodes to highlight branches.\n[Shift+Click] on tree nodes to rotate branches';
			yShift+=14;
		}
		let helpTextY=(colLabelLocation=='top')? paperHeight-yShift : yShift,
			hText=this.paper.text(10,helpTextY,helpText).attr({'text-anchor':'start','font-size':12,fill:'#FFF'}),
			htBox=hText.getBBox(),
			htb=this.paper.rect(htBox.x-6,htBox.y-4,htBox.width+12,htBox.height+8,10).attr({stroke:'none',fill:'#000'}); //,'fill-opacity':0.8
		hText.toFront();
		helpPopup=this.paper.set(htb,hText).attr({'fill-opacity':0,'stroke-opacity':0}).hide().click(showHideHelpText);

		/*** Drawing labels ***/
		function drawLabels() {
//console.time('DRAW_LABELS');
		    /** Rows **/
		    let rowX,rowLabelAnchor;
		    if (rowLabelLocation=='left') {
				rowX=leftSpace+rowLabelSpace-5;
				rowLabelAnchor='end';
		    }
		    else { // right
				rowX=leftSpace+treeSpace.row+mapW+5;
				rowLabelAnchor='start';
		    }
		    if (colLabelLocation=='bottom' && colLabelOrientation=='diagonal') {rowX+=colDiagSpace;}
		    let rowY=topSpace+Math.round(hCellH/2);
		    if (colLabelLocation=='top') {rowY+=colLabelSpace;}
		    else {rowY+=treeSpace.column;} // 0 or 200
			for (let i=0; i<rowList.length; i++) {
				let label;
				if (hCellH >= 5) {
					label=HM.paper.text(rowX,rowY,rowList[i].label)
					.attr({'font-family':rowLabelFontFamily,'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':rowLabelAnchor,'fill':colorList[rowList[i].selectColorIdx]})
					.data('row',rowList[i])
					.hover(function(){setLabelEmphasis(this,'on');},function(){setLabelEmphasis(this,'off');})
					.click(function(e){if (dragContext) {dragContext=null;} else {var l=this; setTimeout(function(){labelClick(l,e);},300);}}) //setTimeout(labelClick,300,this,e) not compatible with IE9
					.dblclick(function(e){labelDoubleClick(this,e);});
					if (rowLabelFontFamily) {label.attr({'font-family':rowLabelFontFamily});}
					if (movableLabel.row) {label.drag(dragLabelMove,dragLabelStart,dragLabelStop);}
				}
				else { // no visible labels
					label=HM.paper.text(rowX,rowY,'').data('row',rowList[i]);
				}
				rowList[i].aspect=label;
				if (i===0) {labelPosRanges.minY=label.attr('y');}
				else if (i==rowList.length-1) {labelPosRanges.maxY=label.attr('y');}
				rowY+=hCellH;
		    }
			/* Row label groups */
			let shX=(rowLabelLocation=='left')? 5 : -5,
				shY=hCellH/2;
			for (let i=0; i<groupList.row.length; i++) {
				let gr=groupList.row[i],
					startY=rowList[gr.labelIndexRange[0]].aspect.attr('y'),
					endY=rowList[gr.labelIndexRange[1]].aspect.attr('y'),
					b=rowList[gr.labelIndexRange[2]].aspect.getBBox(), //longest label
					barX,grX;
				if (rowLabelLocation=='left') {
					barX=b.x - 7.5;
					grX=barX - 7.5;
				}
				else {
					barX=b.x + b.width + 7.5;
					grX=barX + 7.5;
				}
				let bar=HM.paper.path('M'+(barX+shX)+','+(startY-shY)+' L'+barX+','+startY+' L'+barX+','+endY+' L'+(barX+shX)+','+(endY+shY)).attr({stroke:'#000','stroke-width':1}),
					grLabel=HM.paper.text(grX,startY+((endY-startY)/2),gr.label)
					.attr({'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':rowLabelAnchor,'fill':'#000'})
					.data({group:gr,bar:bar})
					.hover(function(){setLabelEmphasis(this,'on');},function(){setLabelEmphasis(this,'off');});
				if (rowGroupLabelFontFamily) {grLabel.attr({'font-family':rowGroupLabelFontFamily});}
				if (rowGroupOnClick) {grLabel.click(function(e){rowGroupOnClick(this.data('group').id);});}
			}

//console.log('minY='+labelPosRanges.minY+', maxY='+labelPosRanges.maxY);
		    /** Columns **/
		    let colX=leftSpace+Math.round(hCellW/2);
		    if (rowLabelLocation=='left') {colX+=rowLabelSpace;}
		    else {colX+=treeSpace.row;} // 0 or 200
		    let colY,colLabelAnchor;
		    if (colLabelLocation=='top') {
				colY=topSpace+colLabelSpace-5;
				colLabelAnchor='start';
		    }
		    else { // bottom
				colY=topSpace+treeSpace.column+mapH+5;
				colLabelAnchor='end';
				if (colLabelOrientation=='diagonal') {
					if (rowLabelLocation=='right') {colX+=colDiagSpace;}
				}
		    }
		    for (let i=0; i<columnList.length; i++) { // columns
				let label=HM.paper.text(colX,colY,columnList[i].label).
				attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':colLabelAnchor,'fill':colorList[columnList[i].selectColorIdx]})
				//.rotate(colLabelAngle,colX,colY)
				.transform('r'+colLabelAngle+','+colX+','+colY)
				.data('column',columnList[i])
				.hover(function(){setLabelEmphasis(this,'on');},function(){setLabelEmphasis(this,'off');})
				.click(function(e){if (dragContext) {dragContext=null;} else {var l=this; setTimeout(function(){labelClick(l,e);},300);}}) //setTimeout(labelClick,300,this,e) not compatible with IE9
				.dblclick(function(e){labelDoubleClick(this,e);});
				if (colLabelFontFamily) {label.attr({'font-family':colLabelFontFamily});}
				if (movableLabel.column) label.drag(dragLabelMove,dragLabelStart,dragLabelStop);
				columnList[i].aspect=label;
//if (i==0) console.log('DRAW: '+columnList[i].aspect.transform(),' <- ','r'+colLabelAngle+','+colX+','+colY);  //.toString()
				if (i===0) {labelPosRanges.minX=label.attr('x');} // rotation x->y
				else if (i==columnList.length-1) {labelPosRanges.maxX=label.attr('x');}
				colX+=hCellW;
		    }
			/* Column label groups */
			shX=hCellW/2; shX2=shY=null;
			if (colLabelLocation=='top') {
				shX2=(colLabelOrientation=='diagonal')? 7.5 : 0;
				shY=5;
			}
			else {
				shX2=(colLabelOrientation=='diagonal')? -7.5 : 0;
				shY=-5;
			}
//console.log(colLabelFontSize,columnList[groupList.column[0].labelIndexRange[0]].aspect.attr('font-size'));
			for (let i=0; i<groupList.column.length; i++) {
				let gr=groupList.column[i],
					b=columnList[gr.labelIndexRange[2]].aspect.getBBox(), //longest label
//HM.paper.rect(b.x,b.y,b.width,b.height);
					startX,endX;
				if (colLabelOrientation=='diagonal') {
					//var bs=columnList[gr.labelIndexRange[0]].aspect.getBBox();
					//var be=columnList[gr.labelIndexRange[1]].aspect.getBBox();
					if (colLabelLocation=='top') { // use longest label and shift left/right to find start/end
						//startX=bs.x2+7.5; // x2 is right in label orientation
						//endX=be.x2+7.5;
						startX=b.x2-(hCellW*(columnList[gr.labelIndexRange[2]].rank-columnList[gr.labelIndexRange[0]].rank))+7.5; // x2 is right in label orientation
						endX=b.x2+(hCellW*(columnList[gr.labelIndexRange[1]].rank-columnList[gr.labelIndexRange[2]].rank))+7.5;
					}
					else {
						//startX=bs.x-7.5; // x is right in label orientation
						//endX=be.x-7.5;
						startX=b.x-(hCellW*(columnList[gr.labelIndexRange[2]].rank-columnList[gr.labelIndexRange[0]].rank))-7.5; // x2 is right in label orientation
						endX=b.x+(hCellW*(columnList[gr.labelIndexRange[1]].rank-columnList[gr.labelIndexRange[2]].rank))-7.5;
					}
				}
				else {
					startX=columnList[gr.labelIndexRange[0]].aspect.attr('x');
					endX=columnList[gr.labelIndexRange[1]].aspect.attr('x');
				}
				let barY,grY;
				if (colLabelLocation=='top') {
					barY=b.y - 7.5;
					grY=barY - 7.5;
				}
				else {
					barY=b.y2 + 7.5;
					grY=barY + 7.5;
				}
				let bar=HM.paper.path('M'+(startX-shX-shX2)+','+(barY+shY)+'L'+startX+','+barY+' L'+endX+','+barY+' L'+(endX+shX-shX2)+','+(barY+shY)).attr({stroke:'#000','stroke-width':1}),
					grX=startX+(shX2)+(endX-startX)/2,
					grLabel=HM.paper.text(grX,grY,gr.label)
					.attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':colLabelAnchor,'fill':'#000'})
					.transform('r'+colLabelAngle+','+grX+','+grY)
					.data({group:gr,bar:bar})
					.hover(function(){setLabelEmphasis(this,'on');},function(){setLabelEmphasis(this,'off');});
				if (colGroupLabelFontFamily) {label.attr({'font-family':colGroupLabelFontFamily});}
				if (columnGroupOnClick) {grLabel.click(function(e){columnGroupOnClick(this.data('group').id);});}
			}

		}
    };


	/*************************************************************/
    /********************** Computing heatmap ********************/
    /*************************************************************/
	function applyMaxRange() {
		function computeQuartile(ordValues,pos) {
			let calcValue;
			if (pos==Math.floor(pos)) {calcValue=ordValues[pos-1];}
			else {
				pos=Math.floor(pos);
				calcValue=(ordValues[pos-1] + ordValues[pos]) / 2;
			}
			return calcValue;
		}
		if (maxRange.type.match(/^abs/)) { // absolute
			rangeMinValue=maxRange.value[0];
			rangeMaxValue=maxRange.value[1];
		}
		else if (maxRange.type.match(/^qua/)) { // quartile
			//Compute quartiles
			let orderedValues=[];
			for (let i=0; i<rowList.length; i++) {
				for (let j=0; j<rowList[i].heatCellList.length; j++) {
					if (rowList[i].heatCellList[j].value !== null) orderedValues.push(rowList[i].heatCellList[j].value);
				}
			}
			orderedValues.sort(idvis.sortNumber);
			var lowQuartValue=computeQuartile(orderedValues,orderedValues.length / 4),
				midQuartValue=computeQuartile(orderedValues,orderedValues.length / 2),
				upQuartValue=computeQuartile(orderedValues,3 * orderedValues.length / 4);
			rangeMinValue=midQuartValue - (midQuartValue - lowQuartValue) * maxRange.value;
			rangeMaxValue=midQuartValue + (upQuartValue - midQuartValue) * maxRange.value;
//console.log(lowQuartValue,midQuartValue,upQuartValue);
		}
//console.log(rangeMinValue,rangeMaxValue);
		// Update cell values
		for (let i=0; i<rowList.length; i++) {
			for (let j=0; j<rowList[i].heatCellList.length; j++) {
				if (rowList[i].heatCellList[j].value !== null) {
					if (rowList[i].heatCellList[j].value < rangeMinValue) {rowList[i].heatCellList[j].value=rangeMinValue;}
					else if (rowList[i].heatCellList[j].value > rangeMaxValue) {rowList[i].heatCellList[j].value=rangeMaxValue;}
				}
			}
		}
	}

    this.computeHeatMap=function() {
//console.time('COMPUTE');
		/* Reset some palette properties (because of possible switch from 2-areas to 1-area palette) */
		palette.hasNegValues=false;
		palette.hasPosValues=false;
		if (palette.trimValues.svgLeft) {palette.trimValues.svgLeft.data('line').remove(); palette.trimValues.svgLeft.remove();}
		if (palette.trimValues.svgRight) {palette.trimValues.svgRight.data('line').remove(); palette.trimValues.svgRight.remove();}
        palette.trimValues={left:0,right:0,svgLeft:null,svgRight:null};

		if (HM.normProcess.reference=='z-score') {
			if (!maxAbsZscore) {
				processData();
//console.log('Max z-score='+maxAbsZscore);
			}
			/* Convert Z-scores to % */
			//convertZscores2pc();
			for (let i=0; i<rowList.length; i++) {
				for (let j=0; j<rowList[i].heatCellList.length; j++) {
					if (rowList[i].heatCellList[j].value !== null && !rowList[i].heatCellList[j].isExcluded) {
						rowList[i].heatCellList[j].pcValue=(rowList[i].heatCellList[j].zscore===null)? 100 : Math.round(100*(rowList[i].heatCellList[j].zscore/maxAbsZscore)); //reference=0
						if (rowList[i].heatCellList[j].pcValue < 0) {palette.hasNegValues=true;}
						else {palette.hasPosValues=true;}
					}
				}
			}
		}
		else {processData();}
        
        if (palette.hasNegValues) palette.trimValues.left=100;
        if (palette.hasPosValues) palette.trimValues.right=100;
        
		this.computeCellsColor();
//console.timeEnd('COMPUTE');
    };

    function processData() {
	    if (HM.normProcess.scope=='row') {
		    for (let i=0; i<rowList.length; i++) {
			    normalizeData(rowList[i].heatCellList);
		    }
	    }
	    else if (HM.normProcess.scope=='column') {
		    for (let i=0; i<columnList.length; i++) {
			    normalizeData(columnList[i].heatCellList);
		    }
	    }
	    else { // entire data set
		    let cellList=[];
		    for (let i=0; i<rowList.length; i++) {
			    for (let j=0; j<rowList[i].heatCellList.length; j++) {
				    cellList.push(rowList[i].heatCellList[j]);
			    }
		    }
		    normalizeData(cellList);
	    }
    }

    function normalizeData(cellList) {
		let valueList=[];
		for (let i=0; i<cellList.length; i++) {
			if (cellList[i].value !== null && !cellList[i].isExcluded) valueList.push(cellList[i].value);
		}
		if (HM.normProcess.reference=='z-score') {
			let mean=0;
			for (let i=0; i<valueList.length; i++) {
				mean+=valueList[i];
			}
			mean/=valueList.length; // mean value
			let sqDelta=0;
			for (let i=0; i<valueList.length; i++) {
				sqDelta+=Math.pow(valueList[i]-mean,2); // delta^2
			}
			let stdDev=Math.sqrt(sqDelta/valueList.length);
//console.log('StdDev='+stdDev);

			for (let i=0; i<cellList.length; i++) { // compute z-scores & record biggest
				if (cellList[i].value !== null && !cellList[i].isExcluded) {
					if (stdDev) {
						cellList[i].zscore=(cellList[i].value-mean)/stdDev; // z-score: (x-mean)/stdDev
						maxAbsZscore=Math.max(maxAbsZscore,Math.abs(cellList[i].zscore));
					}
					else {cellList[i].zscore=null;}
//console.log('maxZ='+maxAbsZscore+' / '+cellList[i].zscore);
			    //valueMaxRange=Math.max(Math.abs(cellList[i].zscore),valueMaxRange);
				}
			}
			/*
			for (var i=0; i<cellList.length; i++) { // z-score -> pcValue
				if (cellList[i].value != null) {
					cellList[i].pcValue=Math.round(100*cellList[i].zscore/valueMaxRange);
				}
			}
			*/
//console.log(HM.normProcess.reference+': Mean='+mean+', Range='+valueMaxRange);
		}
		else {
			valueList.sort(idvis.sortNumber);
			let numValues=valueList.length,
				minValue=valueList[0],
				maxValue=valueList[numValues-1],
//for (var i=0; i<valueList.length; i++) {console.log(valueList[i]);}
				referenceValue=0;
			if (HM.normProcess.reference=='mean') {
				for (let i=0; i<valueList.length; i++) {referenceValue+=valueList[i];}
				referenceValue/=valueList.length;
			}
			if (HM.normProcess.reference=='median') {
				let halfNumValues=numValues/2;
				referenceValue=(Math.floor(halfNumValues)==halfNumValues)? (valueList[halfNumValues]+valueList[halfNumValues-1])/2 : valueList[Math.floor(halfNumValues)]; // indexes!!!
			}
			else if (HM.normProcess.reference=='min') {referenceValue=minValue;}
			else if (HM.normProcess.reference=='zero') {referenceValue=0;}
			else if (HM.normProcess.reference=='mid-range') {referenceValue=minValue+(maxValue-minValue)/2;}
			else if (HM.normProcess.reference=='user') {
				if (typeof(HM.normProcess.limitValues.min) !== 'undefined') minValue=HM.normProcess.limitValues.min;
				if (typeof(HM.normProcess.limitValues.max) !== 'undefined') maxValue=HM.normProcess.limitValues.max;
				referenceValue=(typeof(HM.normProcess.limitValues.ref) !== 'undefined')? HM.normProcess.limitValues.ref : (maxValue-minValue)/2;
			}
			var valueMaxRange=Math.max(referenceValue-minValue,maxValue-referenceValue);
//console.log(referenceValue+','+valueMaxRange);
			for (let i=0; i<cellList.length; i++) {
				if (cellList[i].value !== null) {
					cellList[i].pcValue=Math.round(100*((cellList[i].value-referenceValue)/valueMaxRange));
					if (cellList[i].value < referenceValue) {palette.hasNegValues=true;}
					else {palette.hasPosValues=true;}
				}
			}
		}
    }

//    function convertZscores2pc() {
//		for (var i=0; i<rowList.length; i++) {
//			for (var j=0; j<rowList[i].heatCellList.length; j++) {
//				if (rowList[i].heatCellList[j].value != null && !rowList[i].heatCellList[j].isExcluded) {
//					rowList[i].heatCellList[j].pcValue=(rowList[i].heatCellList[j].zscore==null)? 100 : Math.round(100*(rowList[i].heatCellList[j].zscore/maxAbsZscore)); //reference=0
//				}
//			}
//		}
//    }

    this.computeCellsColor=function(cellList,noDistrib) {
		//cellHighlightColor=(this.normProcess.colors.match('.b.'))? '#FF0' : '#0DD'; // black as middle color (this.normProcess.colors.match('BW|BY'))? '#0DD' : '#0AA';
		cellHighlightColor=(!this.normProcess.colors.match('Y|W'))? '#FF0' : (!this.normProcess.colors.match('R'))? '#E0E' : '#0DD'; // black as middle color (this.normProcess.colors.match('BW|BY'))? '#0DD' : '#0AA';
		noDataCellColor=(this.normProcess.colors.match('W'))? '#000' : '#EEE';
		if (!cellList) {
			cellList=[];
			for (let i=0; i<rowList.length; i++) {
				for (let j=0; j<rowList[i].heatCellList.length; j++) {
					cellList.push(rowList[i].heatCellList[j]);
				}
			}
		}
		let selColors=HM.normProcess.colors.split(''),
			downColorCode=selColors[0]+selColors[1],
			upColorCode=selColors[1]+selColors[2],
            scaleTrimLeft,scaleTrimRight;
        if (palette.hasNegValues && palette.hasPosValues) {
            scaleTrimLeft=100/palette.trimValues.left;
            scaleTrimRight=100/palette.trimValues.right;
        }
        else if (palette.hasNegValues) {
            scaleTrimLeft=100/(palette.trimValues.left-palette.trimValues.right);
        }
        else { // has only pos values
            scaleTrimRight=100/(palette.trimValues.right-palette.trimValues.left);
		}
		// Circle? Compute value2 range
		var cellMaxRadius;
		if (roundedCells) {
			if (!idvis.isNumber(roundedCells.range[0]) || !idvis.isNumber(roundedCells.range[1])) {
				let [minValue2,maxValue2]=[null,null];
				for (let i=0; i<cellList.length; i++) {
					const hcValue=(roundedCells.fitMainValue)? cellList[i].value : cellList[i].value2;
					if (idvis.isNumber(hcValue)) {
						if (minValue2===null) {
							minValue2=maxValue2=hcValue;
						}
						else {
							minValue2=Math.min(minValue2,hcValue);
							maxValue2=Math.max(maxValue2,hcValue);
						}
					}
				}
				if (!idvis.isNumber(roundedCells.range[0])) {roundedCells.range[0]=minValue2;}
				if (!idvis.isNumber(roundedCells.range[1])) {roundedCells.range[1]=maxValue2;}
			}
			cellMaxRadius=Math.sqrt((roundedCells.range[1]-roundedCells.range[0])/Math.PI); // surface -> radius
		}

		for (let i=0; i<cellList.length; i++) {
			let hCell=cellList[i];
			if (hCell.isExcluded) {hCell.color='#999';}
			else if (hCell.value === null) {hCell.color=noDataCellColor;}
			else {
				if (palette.hasNegValues && palette.hasPosValues) {
                    if (hCell.pcValue < 0) { // below reference
                        //let pcVal=Math.min(100,Math.abs(hCell.pcValue)), //Math.round(100*((referenceValue-hCell.value)/valueMaxRange));
                        let usedValueUp=(Math.min(100,-hCell.pcValue*scaleTrimLeft)),
                            usedValueDown=Math.max(0,100+hCell.pcValue*scaleTrimLeft),
                            tempColorStrg=colorPairSet.DOWN[downColorCode].replace(/\-/g,usedValueDown);
                        hCell.color=tempColorStrg.replace(/\+/g,usedValueUp);
                    }
                    else { // above reference
                        let usedValueUp=Math.max(0,(Math.min(100,hCell.pcValue*scaleTrimRight))),
                            usedValueDown=100-usedValueUp,
                            tempColorStrg=colorPairSet.UP[upColorCode].replace(/\-/g,usedValueDown);
                        hCell.color=tempColorStrg.replace(/\+/g,usedValueUp);
                    }
                }
                else if (palette.hasNegValues) { // neg values only
                    let usedValueUp=Math.max(0,(Math.min(100,(-hCell.pcValue-palette.trimValues.right)*scaleTrimLeft))),
                        usedValueDown=100-usedValueUp, 
                        tempColorStrg=colorPairSet.DOWN[downColorCode].replace(/\-/g,usedValueDown);
                    hCell.color=tempColorStrg.replace(/\+/g,usedValueUp);
                }
                else { // Pos values only
                    let usedValueUp=Math.max(0,(Math.min(100,(hCell.pcValue-palette.trimValues.left)*scaleTrimRight))),
                        usedValueDown=100-usedValueUp,
                        tempColorStrg=colorPairSet.UP[upColorCode].replace(/\-/g,usedValueDown);
                    hCell.color=tempColorStrg.replace(/\+/g,usedValueUp);
                }
			}
			mapContext.fillStyle=hCell.color;
			//mapContext.globalAlpha = (hCell.isFlagged)? 0.2 : 1;
			if (roundedCells) { // circle
if (hCell.value !== null) {
				const centerX=Math.round(hCellW*hCell.columnIdx+midCell),
					  centerY=Math.round(hCellH*hCell.rowIdx+midCell),
					  hcValue=(roundedCells.fitMainValue)? hCell.value : hCell.value2;
					  value=Math.min(Math.max(hcValue,roundedCells.range[0]),roundedCells.range[1]), // make sure values are within range
					  radius=(hcValue)? Math.max(1,Math.round(defaultRadius*( Math.sqrt((value-roundedCells.range[0])/Math.PI) / cellMaxRadius))) : defaultRadius;
				mapContext.beginPath();
				mapContext.arc(centerX, centerY, radius, 0, 2 * Math.PI, false);
				mapContext.fill();
}
			}
			else { // default: rectangle
				mapContext.fillRect(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx,hCellW,hCellH);
			}
			//Flag
			if (hCell.isFlagged && cellFlagVisible) {drawCellFlag(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx);}
			//Selection
			if (hCell.isSelectedData) {hCell.isSelectedData.attr({stroke:cellHighlightColor});}
		}

		/* Cell value trace */
		let chkTrace=document.getElementById(formDivID+'_trace');
		if (chkTrace) {
			if ((HM.normProcess.scope=='column' && hCellW < 20) || (HM.normProcess.scope != 'column' && hCellH < 20)) {chkTrace.disabled=true;}
			else {
				chkTrace.disabled=false;
				if (cellsValueTrace) {drawCellsValueTrace(cellList);}
			}
		}

		drawPalette(noDistrib);

		if (hCellFrame) {hCellFrame.attr({stroke:cellHighlightColor});} // update cell highlight frame if already initialized
		if (focusArea) {  // update zoom area if already initialized
			focusArea.attr({stroke:cellHighlightColor}); // ,fill:cellHighlightColor
			focusArea.data('close')[0].attr({fill:cellHighlightColor});
			if (focusArea.data('selRows')) {focusArea.data('selRows')[0].attr({fill:cellHighlightColor});}
			if (focusArea.data('selColumns')) {focusArea.data('selColumns')[0].attr({fill:cellHighlightColor});}
		}
	};

	this.updateCellsValue=function(scope,itemId,updateString) {
		let labelIdxList={row:{},column:{}}, // type.id -> idx
			selLabel;
		//for (let type in labelIdxList) {
        Object.keys(labelIdxList).forEach(function(type) {
			let labelList=(type=='row')? rowList : columnList;
			for (let i=0; i<labelList.length; i++) {
				labelIdxList[type][labelList[i].id]=i;
				if (type==scope && labelList[i].id==itemId) {selLabel=labelList[i];}
			}
		});
		let updateData=updateString.split(';');
		for (let i=0; i<updateData.length; i++) {
			let cellData=updateData[i].split(','),
				rowIdx,colIdx;
			if (scope=='row') {
				rowIdx=labelIdxList.row[itemId];
				colIdx=labelIdxList.column[cellData[0]];
			}
			else { // column
				rowIdx=labelIdxList.row[cellData[0]];
				colIdx=labelIdxList.column[itemId];
			}
			let hCell=rowList[rowIdx].heatCellList[colIdx];
			hCell.value=(cellData[1]===null || cellData[1].length===0 || isNaN(cellData[1]*1))? null : cellData[1]*1; // null*1=''*1=0!!!!
			//if (hCell.value==null && hCell.isSelectedData) {selectCell(hCell);} //unselect cell
		}
		if (scope==HM.normProcess.scope && HM.normProcess.reference != 'z-score') { // update a single row or column
			normalizeData(selLabel.heatCellList);
			HM.computeCellsColor(selLabel.heatCellList);
		}
		else { // all other cases
			maxAbsZscore=0; // reset z-score calculation
			HM.computeHeatMap();
		}
	};

	function updateCellsExclusion(label,excludedStatus) { // label object,true/false
		for (let i=0; i<label.heatCellList.length; i++) {
			label.heatCellList[i].isExcluded=excludedStatus;
		}
		maxAbsZscore=0; // reset z-score calculation
		if (HM.normProcess.scope=='all' || HM.normProcess.reference=='z-score' || label.type != HM.normProcess.scope) {
			HM.computeHeatMap();
		}
		else {
			normalizeData(label.heatCellList);
			HM.computeCellsColor(label.heatCellList);
		}
	}

	/*** Palette ***/
	function drawPalette(noDistrib) {
		/* Compute values distribution */
		if (!noDistrib) {
			let binSize,halfRange;
			if (palette.hasNegValues && palette.hasPosValues) {
				binSize=5;
				halfRange=100/binSize;
			}
			else {
				binSize=2.5;
				halfRange=(palette.hasNegValues)? 100/binSize : 0;
			}
			let binHalf=binSize/2;
			for (let b=0; b<=40; b++) {palette.valueDistrib[b]=0;} // 41 bins
			for (let i=0; i<rowList.length; i++) {
				for (let j=0; j<rowList[i].heatCellList.length; j++) {
					let hCell=rowList[i].heatCellList[j];
					if (hCell.isExcluded || hCell.value===null) {continue;}
					let bin;
					if (hCell.pcValue < 0) {
						let pcValue=Math.min(100,Math.abs(hCell.pcValue));
						bin=halfRange-Math.floor((pcValue+binHalf)/binSize);
					}
					else {
						let pcValue=Math.min(100,hCell.pcValue);
						bin=halfRange+Math.floor((pcValue+binHalf)/binSize);
					}
					palette.valueDistrib[bin]++;
				}
			}
		}
//console.log(palette.valueDistrib);

		/* (Re)Draw palette */
		let palSVG=palette.svgElements,
            palDistr=palette.valueDistrib,
			palX=palette.posX+0.5, palY=palette.posY+0.5;
		//Clear palette
		if (palSVG.length) {
			for (let i=0; i<palSVG.length; i++) {palSVG[i].remove();}
		}
		palSVG.length=0;

		const palH=100, palW=200, palM=palW/2; // palette dimensions
		//Draw color gradients
		let grad=HM.normProcess.colors.split(''), // 3 letter-array
        palTrim=palette.trimValues;
		if (palette.hasNegValues && palette.hasPosValues) {
			//Negative values
			palSVG.push(HM.paper.rect(palX+palM-palTrim.left,palY,palTrim.left,palH).attr({stroke:'none',fill:'180-'+colorCharCode[grad[1]]+'-'+colorCharCode[grad[0]]})); // gradient
			palSVG.push(HM.paper.rect(palX,palY,palM-palTrim.left,palH).attr({stroke:'none',fill:colorCharCode[grad[0]]}));
			//Positive values
			palSVG.push(HM.paper.rect(palX+palM,palY,palTrim.right,palH).attr({stroke:'none',fill:'0-'+colorCharCode[grad[1]]+'-'+colorCharCode[grad[2]]})); // gradient
			palSVG.push(HM.paper.rect(palX+palM-1+palTrim.right,palY,palM-palTrim.right,palH).attr({stroke:'none',fill:colorCharCode[grad[2]]}));
			palette.refLimit=palX+palM;
		}
		else {
			if (palette.hasNegValues) {
                let range=(palTrim.left-palTrim.right)*2;
				palSVG.push(HM.paper.rect(palX-1+palW-palTrim.right*2,palY,palTrim.right*2,palH).attr({stroke:'none',fill:colorCharCode[grad[1]]})); // right
                palSVG.push(HM.paper.rect(palX-1+palW-palTrim.right*2-range+1,palY,range,palH).attr({stroke:'none',fill:'180-'+colorCharCode[grad[1]]+'-'+colorCharCode[grad[0]]})); // middle
				palSVG.push(HM.paper.rect(palX,palY,palW-range-palTrim.right*2,palH).attr({stroke:'none',fill:colorCharCode[grad[0]]})); // left
			}
			else {
				let range=(palTrim.right-palTrim.left)*2;
                palSVG.push(HM.paper.rect(palX,palY,palTrim.left*2,palH).attr({stroke:'none',fill:colorCharCode[grad[1]]})); // left
				palSVG.push(HM.paper.rect(palX-1+palTrim.left*2,palY,range+1,palH).attr({stroke:'none',fill:'0-'+colorCharCode[grad[1]]+'-'+colorCharCode[grad[2]]})); // middle
				palSVG.push(HM.paper.rect(palX-1+palTrim.right*2,palY,palW-palTrim.right*2,palH).attr({stroke:'none',fill:colorCharCode[grad[2]]})); // right
			}
            palette.refLimit=null;
		}
		//Draw histogram of distribution
		let binMaxValue=0;
		for (let i=0; i<palDistr.length; i++) {binMaxValue=Math.max(binMaxValue,palDistr[i]);}
		let distScaleH=palH/(1.1*binMaxValue),
			dpStrg='M'+palX+','+(palY+palH),
			curY=palY;
		for (let i=0; i<=40; i++) { // 41 bins
			let dx=(i===0 || i==40)? 3 : 5,
				y=palY-Math.round(palDistr[i]*distScaleH);
			dpStrg+=' l0,'+(y-curY)+' l'+dx+',0';
			curY=y;
		}
		dpStrg+=' l0,'+(palY-curY);
		palSVG.push(HM.paper.path(dpStrg).attr({stroke:cellHighlightColor,'stroke-width':1}));
		//Draw reference (if any)
		if (palette.hasNegValues && palette.hasPosValues) {
			palSVG.push(HM.paper.path('M'+(palX+100)+','+palY+' l0,'+palH).attr({stroke:cellHighlightColor,'stroke-width':1}));
		}
		//Draw frame
		palSVG.push(HM.paper.rect(palX,palY,palW,palH));
		//Trim handles
        let scale=(palette.hasNegValues && palette.hasPosValues)? 1 : 2;
		let trim=palette.trimValues;
        /* Min */
        if (trim.svgLeft) {trim.svgLeft.data('line').attr({stroke:cellHighlightColor});}
        else {
            let tX=(palette.refLimit)? palX+(100-trim.left)*scale : (palette.hasPosValues)? palX+trim.left*scale : palX+(100-trim.left)*scale,
                line=HM.paper.path('M'+tX+','+palY+' l0,100').attr({stroke:cellHighlightColor,'stroke-width':2}).hide();
            trim.svgLeft=HM.paper.path('M'+tX+','+(palY+101)+' l8,8 l-17,0 Z').attr({fill:'#000',stroke:'none'})
            .data({type:'left',pos:tX,totShift:0,shift:0,line:line}).hover(function(){this.attr({fill:'#F00'});},function(){if(!dragContext){this.attr({fill:'#000'});}})
            .drag(dragTrimMove,dragTrimStart,dragTrimStop);
        }
        /* Max */
        if (trim.svgRight) {trim.svgRight.data('line').attr({stroke:cellHighlightColor});}
        else {
            let tX=(palette.refLimit)? palX+palW-(100-trim.right)*scale : (palette.hasPosValues)? palX+trim.right*scale : palX+palW-trim.right*scale,
                line=HM.paper.path('M'+tX+','+palY+' l0,100').attr({stroke:cellHighlightColor,'stroke-width':2}).hide();
            trim.svgRight=HM.paper.path('M'+tX+','+(palY+101)+' l8,8 l-17,0 Z').attr({fill:'#000',stroke:'none'})
            .data({type:'right',pos:tX,totShift:0,shift:0,line:line}).hover(function(){this.attr({fill:'#F00'});},function(){if(!dragContext){this.attr({fill:'#000'});}})
            .drag(dragTrimMove,dragTrimStart,dragTrimStop);
        }
    }


    /*********************************************************/
    /********************** Cell events **********************/
    /*********************************************************/
	function findHeatCell(e) {
		//var mouseXY=getMousePositionInChart(e); //chartLibrary2.js function
		//var mouseX,mouseY; // pos in hCellsSensor
		//if (e.offsetX) { // Firefox & IE 9+ (no need to check for IE < 9 because not compatible with CANVAS)
		//	mouseX=e.offsetX; //-cellX0;
		//	mouseY=e.offsetY; //-cellY0;
		//}
		//else {
		//	mouseX=e.layerX-cellX0;
		//	mouseY=e.layerY-cellY0;
		//}
		let mousePos=idvis.getMousePositionInElement(e), // relative to hCellsSensor NOT chart (not clear why; due to multiple DIVs?)
			colIdx=Math.round((mousePos[0]/hCellW)-0.5), // - because index not rank
			rowIdx=Math.round((mousePos[1]/hCellH)-0.5);
//document.getElementById('mX').value=rowIdx;
//document.getElementById('mY').value=colIdx;
		if (rowIdx>=rowList.length) {rowIdx=rowList.length-1;}
		else if (rowIdx<0) {rowIdx=0;}
		if (colIdx>=columnList.length) {colIdx=columnList.length-1;}
		else if (colIdx<0) {colIdx=0;}
		return rowList[rowIdx].heatCellList[colIdx];
	}

    function setCellEmphasis(hCell,action,e) {
		let cellBorderColor,cellCursor, //labelColor,
		//var c=hCell.aspect;
			cX=hCellW*hCell.columnIdx, cY=hCellH*hCell.rowIdx,
			rowLabel=rowList[hCell.rowIdx],colLabel=columnList[hCell.columnIdx];
		if (action=='on') {
			/*
			if (cellPopupText) {
				cellPopupText.remove();
				cellPopupText=null;
			}
			*/
			//labelColor='#F00';
			rowLabel.aspect.attr({fill:'#F00'});
			colLabel.aspect.attr({fill:'#F00'});
			cellBorderColor=cellHighlightColor;
			cellCursor=(cellOnClick || e.ctrlKey || e.shiftKey)? 'pointer' : 'auto';
			let valueStrg;
			if (cellOnMouseOver) {valueStrg=cellOnMouseOver(hCell.value,hCell,rangeMinValue,rangeMaxValue);} // hCell in case more complex display is needed
			else {
				if (hCell.value===null) {valueStrg='* no value *';}
				else {
					let overStrg=(hCell.value <= rangeMinValue)? '≤' : (hCell.value >= rangeMaxValue)? '≥' : '';
					valueStrg='Value: ' + overStrg + idvis.formatNumber(hCell.value,3); //(Math.round(100*hCell.value)/100);
					if (HM.normProcess.reference=='z-score' && !hCell.isExcluded) {
						if (hCell.zscore===null) {valueStrg+=' [no z-score]';}
						else {valueStrg+=' [z='+(Math.round(100*hCell.zscore)/100)+']';}
					}
				}
				if (roundedCells && !roundedCells.fitMainValue) {
					valueStrg+='\n'+roundedCells.sizeLabel+': ';
					valueStrg+=(idvis.isNumber(hCell.value2))? idvis.formatNumber(hCell.value2,3) : '* no value *';
				}
//valueStrg+=' ('+hCell.pcValue+')';
			}
			let textStrg=(e.ctrlKey || e.shiftKey)? HM.entities.row+': '+rowLabel.label+'\n'+HM.entities.column+': '+colLabel.label+'\n'+valueStrg : valueStrg;
			if (hCell.isExcluded) textStrg+='\n* excluded *';
			if (hCell.isFlagged) textStrg+='\n* '+flagText+' *';
			//var textStrg=HM.entities.row+': '+rowLabel.label+'\n'+HM.entities.column+': '+colLabel.label+'\n'+valueStrg;
//textStrg+='\n'+hCell.rowIdx+','+hCell.columnIdx;
			let x=cellX0+cX, // c.attr('x')+Math.round(c.attr('width')/2);
				xm=x+Math.round(hCellW/2),
				y=cellY0+cY; //c.attr('y');
			//drawCellPopup(xm,y,textStrg,cellHighlightColor,cellHighlightBgColor,cellHighlightBgColor);
			drawCellPopup(xm,y,textStrg,'#FF0','#000',cellHighlightColor);
			hCellFrame.attr({x:x+1,y:y+1}).show();
		}
		else { //off
			//labelColor='#000';
			if (rowLabel.selectColorIdx) {rowLabel.aspect.attr({fill:colorList[rowLabel.selectColorIdx]});}
			else {rowLabel.aspect.attr({fill:'#000'});}
			if (colLabel.selectColorIdx) {colLabel.aspect.attr({fill:colorList[colLabel.selectColorIdx]});}
			else {colLabel.aspect.attr({fill:'#000'});}
			cellBorderColor=(hCell.isExcluded)? '#999' : hCell.color; //c.attr('fill');
			cellCursor='auto';
			cellPopupText.remove();
			hCellFrame.hide();
		}
		//mapContext.strokeStyle=cellBorderColor;
		//mapContext.lineWidth=2;
		//mapContext.strokeRect(cX+1,cY+1,hCellW-2,hCellH-2);
		if (hCell.isFlagged && action=='off' && cellFlagVisible) {drawCellFlag(cX,cY);}
    }

	HM.showCellsValueLine=function(status) {
		cellsValueTrace=status;
		let labelList=(HM.normProcess.scope=='column')? columnList : rowList;
		for (let i=0; i<labelList.length; i++) {
			for (let j=0; j<labelList[i].heatCellList.length; j++) {
				let hCell=labelList[i].heatCellList[j];
				mapContext.fillStyle=hCell.color;
				mapContext.fillRect(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx,hCellW,hCellH);
			}
			if (cellsValueTrace) {drawCellsValueTrace(labelList[i].heatCellList);}
		}
		if (existFlaggedCells && cellFlagVisible) {HM.showCellFlag(cellFlagVisible);}
	};

	function drawCellsValueTrace(cellList) {
		if (cellsValueTrace===false) return; // just to be safe
		if ((HM.normProcess.scope=='column' && hCellW < 20) || (HM.normProcess.scope != 'column' && hCellH < 20)) {return;} // cells are too small
		let prevLinePos=null,prevRowIdx=-2,prevColumnIdx=-2, // Not necessarily index in cellList!
			refPos,refScale;
		if (HM.normProcess.scope=='column') {
			if (palette.hasNegValues) {refPos=hCellW/2; refScale=hCellW/2;}
			else {refPos=0; refScale=hCellW;}
		}
		else {
			if (palette.hasNegValues) {refPos=hCellH/2; refScale=hCellH/2;}
			else {refPos=hCellH; refScale=hCellH;}
		}
		for (let i=0; i<cellList.length; i++) {
			let hCell=cellList[i];
			if (hCell.pcValue===null) {
				prevRowIdx=-2;
				prevColumnIdx=-2;
				continue;
			}
			mapContext.strokeStyle=cellHighlightColor;
			mapContext.lineWidth=1;
			mapContext.beginPath();

			let x0=hCellW*hCell.columnIdx,
				y0=hCellH*hCell.rowIdx,
				linePos;
			if (HM.normProcess.scope=='column') {
				linePos=Math.floor(x0+refPos+((hCell.pcValue/100)*refScale))+0.5;
				if (hCell.rowIdx===0) {prevLinePos=null;} // 1st cell in row
				else if (prevColumnIdx != hCell.columnIdx || prevRowIdx != hCell.rowIdx-1) { // recompute prevLinePos
					let prevCell=rowList[hCell.rowIdx-1].heatCellList[hCell.columnIdx]; // get cell in previous row
					prevLinePos=(prevCell.pcValue===null)? null : Math.floor(x0+refPos+((prevCell.pcValue/100)*refScale))+0.5;
				}
				if (prevLinePos===null) {mapContext.moveTo(linePos,y0);}  // No previous cell or it has no value
				else {
					mapContext.moveTo(prevLinePos,y0); //!!!!!!!
					mapContext.lineTo(linePos,y0);
				}
				mapContext.lineTo(linePos,y0+hCellH+1);
//if (i < 10) {console.log(i,linePos);}
			}
			else { // row or both
				linePos=Math.floor(y0+refPos-((hCell.pcValue/100)*refScale))+0.5;
				if (hCell.columnIdx===0) {prevLinePos=null;} // 1st cell in row
				else if (prevRowIdx != hCell.rowIdx || prevColumnIdx != hCell.columnIdx-1) { // recompute prevLinePos
					let prevCell=columnList[hCell.columnIdx-1].heatCellList[hCell.rowIdx]; // get cell in previous column
					prevLinePos=(prevCell.pcValue===null)? null : Math.floor(y0+refPos-((prevCell.pcValue/100)*refScale))+0.5;
				}
				if (prevLinePos===null) {mapContext.moveTo(x0,linePos);}  // No previous cell or it has no value
				else {
					mapContext.moveTo(x0,prevLinePos);
					mapContext.lineTo(x0,linePos);
				}
				mapContext.lineTo(x0+hCellW+1,linePos);
			}

			mapContext.stroke();
			prevRowIdx=hCell.rowIdx;
			prevColumnIdx=hCell.columnIdx;
			prevLinePos=linePos;
		}
/*
		if (HM.normProcess.scope=='column') {
//labelList=,otherDimList,indexType,cellSize
		}
		else { // scope is row or all
			for (let i=0; i<cellList.length; i++) {
				let hCell=cellList[i],
					x0=Math.round(hCellW*hCell.columnIdx),
					y0=hCellH*hCell.rowIdx;
				if (hCell.pcValue===null) {
					prevCellIndex=-2;
					continue;
				}
				mapContext.strokeStyle=cellHighlightColor;
				mapContext.lineWidth=1;
				mapContext.beginPath();
				let linePos=Math.floor(y0+hCellH-((hCell.pcValue/100)*hCellH))+0.5;
				if (prevCellIndex != hCell.columnIdx-1) { // recompute prevLinePos
					if (hCell.columnIdx===0) {prevLinePos=null;} // 1st cell in row
					else {
						let prevCell=columnList[hCell.columnIdx-1].heatCellList[hCell.rowIdx]; // get cell in previous colonne
						prevLinePos=(prevCell.pcValue===null)? null : Math.floor(y0+hCellH-((prevCell.pcValue/100)*hCellH))+0.5;
					}
				}
				if (prevLinePos===null) {mapContext.moveTo(x0,linePos);}  // No previous cell or it has no value
				else {
					mapContext.moveTo(x0,prevLinePos);
					mapContext.lineTo(x0,linePos);
				}
				mapContext.lineTo(x0+hCellW+1,linePos);
				mapContext.stroke();
				prevCellIndex=hCell.columnIdx;
				prevLinePos=linePos;
			}
		}
*/
	}


	function drawCellFlag(cX,cY,color) {
		if (!cellFlagVisible) return; // just to be safe
		if (!color) {color=cellHighlightColor;}
		if (hCellH >= 5) { // top left-side triangle
			mapContext.beginPath();
			mapContext.moveTo(cX,cY);
			mapContext.lineTo(cX+flagSize,cY);
			mapContext.lineTo(cX,cY+flagSize);
			mapContext.lineTo(cX,cY+flagSize);
			mapContext.lineTo(cX,cY);
			mapContext.closePath();
			mapContext.fillStyle=color;
			mapContext.fill();
		}
		else { // left-side rectangle
			mapContext.fillStyle=color;
			mapContext.fillRect(cX,cY,flagSize,hCellH);
		}
	}

	function drawCellPopup(x,y,textStrg,textColor,backColor,borderColor) {
		let t=HM.paper.text(x,y,textStrg).attr({'font-size':10,'text-anchor':'start',fill:textColor}),
			tBox0=t.getBBox();
		t.attr({x:x-Math.round(tBox0.width/2),y:y-Math.round(tBox0.height/2)-5});
		let tBox=t.getBBox(),
			tx=tBox.x,
			ty=tBox.y,
			tw=tBox.width,
			th=tBox.height,
			tb=HM.paper.rect(tx-2,ty-1,tw+4,th+2,4).attr({stroke:borderColor,fill:backColor,'fill-opacity':0.75});
		t.toFront();
		cellPopupText=HM.paper.set(tb,t);
	}

    function selectCell(hCell,forcedStatus,newSelStatus) {
		if (forcedStatus) { // set a specific selection status (independently of cell's current one)
			if (newSelStatus && !hCell.isSelectedData) {
				drawCellSelection(hCell);
			}
			else if (!newSelStatus && hCell.isSelectedData) {
				hCell.isSelectedData.remove();
				hCell.isSelectedData=null;
			}
		}
		else {
			let selStatus;
			if (hCell.isSelectedData) { // unselect
				hCell.isSelectedData.remove();
				hCell.isSelectedData=null;
				selStatus=false;
				if (oneCellSelection) {lastSelectedCell=null;}
			}
			else { // select
				if (oneCellSelection && lastSelectedCell) {selectCell(lastSelectedCell);} // unselect last selected
				drawCellSelection(hCell);
				selStatus=true;
				lastSelectedCell=hCell;
			}
			let midX=cellX0 + (hCell.columnIdx * hCellW) + Math.round(hCellW / 2),
				bottomY=cellY0 + (hCell.rowIdx * hCellH) + hCellH;
			cellOnClick(rowList[hCell.rowIdx].id,columnList[hCell.columnIdx].id,selStatus,midX,bottomY); // user callback (rowId,colId,selStatus,midX,bottomY)
		}
	}

	function drawCellSelection(hCell,shiftX,shiftY) {
		if (!shiftX) shiftX=0;
		if (!shiftY) shiftY=0;
		let s,
			cX=cellX0 + (hCell.columnIdx * hCellW) + Math.round(hCellW / 2) + shiftX,
			cY=cellY0 + (hCell.rowIdx * hCellH) + Math.round(hCellH / 2) + shiftY;
		if (hCellW >=20 && hCellH >=20) {
			let x= cX - 3, y= cY - 2;
			s=HM.paper.path('M'+x+','+y+'l3,5l5,-10').attr({stroke:cellHighlightColor,'stroke-width':3});
		}
		else {
			s=HM.paper.circle(cX,cY,2).attr({stroke:'none',fill:cellHighlightColor});
		}
		hCellsSensor.toFront();
		hCell.isSelectedData=s;
		focusArea.data('close').toFront();
		if (focusArea.data('selRows')) {focusArea.data('selRows').toFront();}
		if (focusArea.data('selColumns')) {focusArea.data('selColumns').toFront();}
    }

    function excludeCell(hCell) {
		maxAbsZscore=0; // reset z-score calculation
		if (hCell.isExcluded) {hCell.isExcluded=false;}
		else {
			hCell.isExcluded=true;
			/*
			if (hCell.isSelectedData) {
				hCell.isSelectedData.remove();
				hCell.isSelectedData=null;
			}
			*/
		}
		if (HM.normProcess.scope=='all' || HM.normProcess.reference=='z-score') {
			HM.computeHeatMap();
		}
		else {
			let label=(HM.normProcess.scope=='row')? rowList[hCell.rowIdx] : columnList[hCell.columnIdx];
			normalizeData(label.heatCellList);
			HM.computeCellsColor(label.heatCellList);
		}
		//HM.colorCells();
		if (cellOnModClick) {
			let midX=cellX0 + (hCell.columnIdx * hCellW) + Math.round(hCellW / 2),
				bottomY=cellY0 + (hCell.rowIdx * hCellH) + hCellH;
			cellOnModClick(rowList[hCell.rowIdx].id,columnList[hCell.columnIdx].id,hCell.isExcluded,midX,bottomY); // user callback (rowId,colId,selStatus,midX,bottomY)
		}
    }

	this.clearAllCells=function() {
		for (let i=0; i<rowList.length; i++) {
			for (let j=0; j<rowList[i].heatCellList.length; j++) {
				let hCell=rowList[i].heatCellList[j];
				if (hCell.isSelectedData) {
					hCell.isSelectedData.remove();
					hCell.isSelectedData=null;
				}
			}
		}
		if (oneCellSelection) {lastSelectedCell=null;}
	};

	this.flagCell=function(colLabelID,rowLabelID,status) {
		//Finding cell coordinates
		let colIdx=-1, rowIdx=-1;
		for (let i=0; i<columnList.length; i++) {
			if (columnList[i].id==colLabelID) {
				colIdx=i;
				break;
			}
		}
		if (colIdx==-1) {return -1;}
		for (let i=0; i<rowList.length; i++) {
			if (rowList[i].id==rowLabelID) {
				rowIdx=i;
				break;
			}
		}
		if (rowIdx==-1) {return -2;}
		let hCell=rowList[rowIdx].heatCellList[colIdx];
		hCell.isFlagged=status;
		let color=(status===true)? cellHighlightColor : (hCell.isExcluded)? '#999' : hCell.color;
		if (cellFlagVisible) drawCellFlag(hCellW*colIdx,hCellH*rowIdx,color);
		return 0;
	};

	this.showCellFlag=function(visStatus) {
		cellFlagVisible=visStatus;
		let flaggedCells=[];
		for (let i=0; i<rowList.length; i++) {
			for (let j=0; j<rowList[i].heatCellList.length; j++) {
				if (rowList[i].heatCellList[j].isFlagged) {
					let hCell=rowList[i].heatCellList[j];
					//canvas
					if (visStatus===false) { // redraw cell w/o flag
						mapContext.fillStyle=hCell.color;
						mapContext.fillRect(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx,hCellW,hCellH);
						flaggedCells.push(hCell);
					}
					else { // show Flag
						drawCellFlag(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx);
					}
					//svg
					if (hCell.isSelectedData) {
						hCell.isSelectedData.remove();
						drawCellSelection(hCell);
					}
				}
			}
		}
		/* Redraw un-flaged cells value line */
		if (cellsValueTrace && flaggedCells.length) {
			drawCellsValueTrace(flaggedCells);
		}
	};

	/**********************************************************/
    /********************* Palette events *********************/
    /**********************************************************/
	function dragTrimStart() {
		dragContext=true;
		this.data('line').toFront().show();
	}
	function dragTrimMove(dx) { // dy is not needed
        let midLimit=(palette.refLimit)? palette.refLimit
                    : (this.data('type')=='left')? palette.trimValues.svgRight.data('pos')+palette.trimValues.svgRight.data('totShift')
                    : palette.trimValues.svgLeft.data('pos')+palette.trimValues.svgLeft.data('totShift'),
//		let midLimit=palette.refLimit,
			shift=this.data('totShift')+dx,
			tX=this.data('pos')+shift;
//console.log(this.data('type'),midLimit,tX);
		if (this.data('type')=='left') {
			if (tX < palette.posX) {tX=palette.posX;}
			else if (tX > midLimit-10) {tX=midLimit-10;}
			shift=tX-this.data('pos');
		}
		else {
			if (tX < midLimit+10) {tX=midLimit+10;}
			else if (tX > palette.posX+200) {tX=palette.posX+200;}
			shift=tX-this.data('pos');
		}
		this.transform('t'+shift+',0');
		this.data('line').transform('t'+shift+',0');
		this.data({shift:shift});
	}
	function dragTrimStop() {
		dragContext=false;
		this.attr({fill:'#000'});
		this.data({totShift:this.data('shift')});
		this.data({Shift:0});
		this.data('line').hide();
		/* Compute new trim value */
		let scale=(palette.hasNegValues && palette.hasPosValues)? 1 : 2;
		if (this.data('type')=='left') {
			palette.trimValues.left=(palette.refLimit)? Math.round((palette.refLimit-this.data('pos')-this.data('totShift'))/scale) : (palette.hasNegValues)? Math.round(100-this.data('totShift')/scale) : Math.round(this.data('totShift')/scale);
            palette.trimValues.left=Math.max(0,Math.min(100,palette.trimValues.left));
        }
		else {
			palette.trimValues.right=(palette.refLimit)? Math.round((this.data('pos')+this.data('totShift')-palette.refLimit)/scale) : (palette.hasPosValues)? Math.round(100+this.data('totShift')/scale) : Math.round(-this.data('totShift')/scale);
            palette.trimValues.right=Math.max(0,Math.min(100,palette.trimValues.right));
        }
		HM.computeCellsColor(null,true); // true to skip computation of values distribution
	}

    /**********************************************************/
    /********************** Label events **********************/
    /**********************************************************/
	function labelClick(l,e) {
		if (dblClickContext) {
//console.log('ABORT CLICK');
			return;
		}
//console.log('*CLICK '+dblClickContext);
		let type=(l.data('row'))? 'row' : 'column',
			returnAction;
		if (type=='row') {
			if (e.ctrlKey || e.shiftKey) {returnAction=(rowOnModClick)? rowOnModClick(l.data(type).id) : null;}
			else {returnAction=(rowOnClick)? rowOnClick(l.data(type).id) : null;}
		}
		else { // column
			if (e.ctrlKey || e.shiftKey) {returnAction=(columnOnModClick)? columnOnModClick(l.data(type).id) : null;}
			else {returnAction=(columnOnClick)? columnOnClick(l.data(type).id) : null;}
		}
		//console.log(returnAction);
		if (returnAction && returnAction.match(/(ex|in)clude/)) {
			let excludedStatus=(returnAction=='exclude')? true : false;
			updateCellsExclusion(l.data(type),excludedStatus);
		}
	}
	function labelDoubleClick(l,e) {
		dblClickContext=true;
//console.log('*DBL-CLICK '+dblClickContext);
		let type=(l.data('row'))? 'row' :'column',
			itemLabel=l.data(type),
			lBox=l.getBBox();
		if (type=='row') {
			let posX=(rowLabelLocation=='left')? lBox.x2 : lBox.x;
			if (e.ctrlKey || e.shiftKey) {if (rowOnModDblClick) {rowOnModDblClick(itemLabel.id,posX,lBox.y);}}
			else if (rowOnDblClick) {rowOnDblClick(itemLabel.id,posX,lBox.y);}
		}
		else { // column
			let posY=(colLabelLocation=='top')? lBox.y2 : lBox.y;
			if (e.ctrlKey || e.shiftKey) {if (columnOnModDblClick) {columnOnModDblClick(itemLabel.id,lBox.x2,posY);}}
			else if (columnOnDblClick) {columnOnDblClick(itemLabel.id,lBox.x2,posY);}
		}
		// Select all cells in row or column
		if (cellOnClick && !oneCellSelection) {
			let newSelStatus=(e.ctrlKey || e.shiftKey)? false : true;
			for (let i=0; i<itemLabel.heatCellList.length; i++) {
				let hCell=itemLabel.heatCellList[i];
				//if ((hCell.isExcluded || hCell.value==null) && newSelStatus) continue; // excluded/no value cells not selected but unselected if previously selected
				if (hCell.isExcluded && newSelStatus) continue; // excluded/no value cells not selected but unselected if previously selected
				selectCell(hCell,true,newSelStatus);
			}
		}
		setTimeout(function(){dblClickContext=false;},500);
	}
    function highLightLabel(label,newColorIdx) {
		//var labelColor=(newColorIdx)? '#00F' : '#000';
		//var labelColor=colorList[newColorIdx];
		label.aspect.attr({fill:colorList[newColorIdx]});
		label.selectColorIdx=newColorIdx;
    }
    function setLabelEmphasis(l,action) {
		let itemLabel,type,isLabel=false,isAnnot=false,isGroup=false,typePopup;
		if (l.data('set')) { // annotation label
			itemLabel=l.data('set');
			type=itemLabel.type;
			isAnnot=true;
			typePopup=(type =='row')? 'column' : 'row';
		}
		else if (l.data('group')) { // group label
			itemLabel=l.data('group');
			type=itemLabel.type;
			isGroup=true;
			typePopup=type;
		}
		else {
			type=(l.data('row'))? 'row' :'column';
			itemLabel=l.data(type);
			isLabel=true;
			typePopup=type;
		}
		let labelColor,labelCursor;
		if (action=='on') {
			if (labelPopupText) return;
			labelColor='#F00';
			labelCursor=(!isLabel)? 'pointer' : (movableLabel[type])? 'move' : 'auto';
			let popupText=itemLabel.popupText || itemLabel.label;
			//if (itemLabel.popupText) {
				if (typePopup=='row') {
					let x=(rowLabelLocation=='left')? l.attr('x')+6 : l.attr('x')-6, //mapGeometry[rowLabelLocation];
						y=l.attr('y'),
					//drawPopup(x,y,itemLabel.popup,'#000','#FFF');
						t=l.paper.text(x,y,popupText).attr({'font-size':12,fill:'#000','text-anchor':'start'}),
						tBox=t.getBBox();
					if (rowLabelLocation=='right') {t.attr('x',t.attr('x')-tBox.width);}
					if (tBox.y < 0) {t.attr('y',t.attr('y')-tBox.y+4);} // tBox.y is negative
					else if (tBox.y+tBox.height > paperHeight) {t.attr('y',t.attr('y')-((tBox.y+tBox.height)-paperHeight+4));}
					tBox=t.getBBox();
					let tx=tBox.x,
						ty=tBox.y,
						tw=tBox.width,
						th=tBox.height,
						tb=l.paper.rect(tx-2,ty-1,tw+4,th+2,4).attr({stroke:'#000',fill:'#FFF','fill-opacity':1});
					t.toFront();
					labelPopupText=l.paper.set(tb,t);
					labelPopupText.transform(l.transform()); // in case inside focusArea
				}
				else { // column
					let x=l.attr('x'),
						y=(colLabelLocation=='top')? l.attr('y')+12 : l.attr('y')-12, //mapGeometry[rowLabelLocation];
						t=l.paper.text(x,y,popupText).attr({'font-size':12,fill:'#000','text-anchor':'middle'}),
						tBox=t.getBBox();
					if (tBox.x < 0) {t.attr('x',t.attr('x')-tBox.x+4);} // tBox.x is negative
					else if (tBox.x+tBox.width > paperWidth) {t.attr('x',t.attr('x')-((tBox.x+tBox.width)-paperWidth+4));}
					tBox=t.getBBox();
					let tx=tBox.x,
						ty=tBox.y,
						tw=tBox.width,
						th=tBox.height,
						tb=l.paper.rect(tx-2,ty-1,tw+4,th+2,4).attr({stroke:'#000',fill:'#FFF','fill-opacity':1});
					t.toFront();
					labelPopupText=l.paper.set(tb,t);
					if (isLabel && hasCells && focusArea.data('status')=='on' && itemLabel.rank > focusArea.data('labelIdx')[2] && itemLabel.rank-1 <= focusArea.data('labelIdx')[3]) { // in case inside focusArea
						if (colLabelLocation=='top') {labelPopupText.transform('t0,'+(focusArea.attr('y')-cellY0));}
						else {labelPopupText.transform('t0,-'+(cellY0+mapH-(focusArea.attr('y')+focusArea.attr('height'))));}
					}
				}
			//}
			/* Group label */
			if (isGroup) { //highlight matching labels
				l.data('bar').attr('stroke',labelColor);
				let labelList=(type=='row')? rowList : columnList;
				for (let i=itemLabel.labelIndexRange[0]; i<=itemLabel.labelIndexRange[1]; i++) {
					labelList[i].aspect.attr({fill:labelColor});
				}
			}
		}
		else { //off
			//labelColor=(itemLabel.selectColorIdx)? '#00F' : '#000';
			labelColor=(isLabel)? colorList[itemLabel.selectColorIdx] : '#000';
			labelCursor='auto';
			if (labelPopupText) {labelPopupText.remove(); labelPopupText=null;}
			/* Group label */
			if (isGroup) { // remove highlight of matching labels
				l.data('bar').attr('stroke','#000');
				let labelList=(type=='row')? rowList : columnList;
				for (let i=itemLabel.labelIndexRange[0]; i<=itemLabel.labelIndexRange[1]; i++) {
					labelList[i].aspect.attr({fill:colorList[labelList[i].selectColorIdx]});
				}
			}
		}
		l.attr({fill:labelColor,cursor:labelCursor});
    }

	/* Multi-labels external action */
	function focussedLabelsAction(call,type) {
		let labelList,selLabelList;
		if (call=='tree') {
			labelList=(type=='row')? columnList : rowList; //row/col inverted for tree
			selLabelList=selectedLabels;
			if (labelsOnSelect[type].ranking) {
				let rankedData=[];
				for (let i=0; i<labelList.length; i++) {rankedData.push(labelList[i].id);}
				labelsOnSelect[type].action(selLabelList,rankedData);
			}
			else {labelsOnSelect[type].action(selLabelList);}
		}
		else { // call=focus
			let targetList,selIdx,idx;
			if (type=='row') {
				targetList=rowList;
				selIdx=0;
				idx=2;
				labelList=columnList;
			}
			else {
				targetList=columnList;
				selIdx=2;
				idx=0;
				labelList=rowList;
			}
			selLabelList=[];
			for (let i=focusArea.data('labelIdx')[selIdx]; i<=focusArea.data('labelIdx')[selIdx+1]; i++) {selLabelList.push(targetList[i].id);}
			if (labelsOnSelect[type].ranking) {
				let rankedData=[];
				for (let i=focusArea.data('labelIdx')[idx]; i<=focusArea.data('labelIdx')[idx+1]; i++) {rankedData.push(labelList[i].id);}
				labelsOnSelect[type].action(selLabelList,rankedData);
			}
			else {labelsOnSelect[type].action(selLabelList);}
		}
	}

    /*** Dragging labels ***/
    function dragLabelStart() { // also fired on click
		dragContext=false;
	}
	function dragLabelInit(label) { // called once by dragLabelMove()
		if (labelPopupText) {labelPopupText.remove(); labelPopupText=null;} //will be recreated by hover during drag
		//dragContext=false;
		dragParams.labelList={};
		if (label.data('row')) {
			dragParams.itemLabel=label.data('row');
			for (let i=0; i<rowList.length; i++) {
				dragParams.labelList[rowList[i].rank]=rowList[i];
			}
			dragParams.maxRank=rowList.length;
			dragParams.type='row';
			dragParams.axis='y';
			dragParams.emptyPosL=label.attr('y');
		}
		else {
			dragParams.itemLabel=label.data('column'); // row or column label
			for (let i=0; i<columnList.length; i++) {
				dragParams.labelList[columnList[i].rank]=columnList[i];
			}
			dragParams.maxRank=columnList.length;
			dragParams.type='column';
			dragParams.axis='x';
			dragParams.emptyPosL=label.attr('x');
		}
		dragParams.movedLabelRanks={min:dragParams.itemLabel.rank,max:dragParams.itemLabel.rank};
		dragParams.startPosL=dragParams.emptyPosL; // original pos at drag start
		dragParams.dragPrevLength=0;
		label.toFront();
		label.attr('fill','#F00');
		dragParams.svgCells=[];

		if (dragParams.itemLabel.heatCellList.length <= 10) { // few cells to move
			for (let i=0; i<dragParams.itemLabel.heatCellList.length; i++) {
				let hCell=dragParams.itemLabel.heatCellList[i],
					svgCell=HM.paper.rect(cellX0+hCellW*hCell.columnIdx+1,cellY0+hCellH*hCell.rowIdx+1,hCellW-2,hCellH-2,0).attr({fill:hCell.color,'fill-opacity':0.5,'stroke-width':2,stroke:cellHighlightColor});
				dragParams.svgCells.push(svgCell);
				if (hCell.isSelectedData) {hCell.isSelectedData.toFront();}
			}
		}
		else { // draw a single svg rect for all cells
			let setX,setY,setW,setH; // row/column of cells selected for drag
			if (dragParams.type=='row') {
				setX=cellX0;
				setY=cellY0+hCellH*(dragParams.itemLabel.rank-1);
				setW=hCellW*columnList.length; //dragParams.itemLabel.heatCellList.length;
				setH=hCellH;
			}
			else { // column
				setX=cellX0+hCellW*(dragParams.itemLabel.rank-1);
				setY=cellY0;
				setW=hCellW;
				setH=hCellH*rowList.length; // dragParams.itemLabel.heatCellList.length;
			}
			dragParams.svgCells.push(HM.paper.rect(setX,setY,setW,setH,0).attr({fill:cellHighlightColor,'fill-opacity':0.5,stroke:'none'}));
		}
//		hCellsSensor.toFront();
//		focusArea.data('close').toFront();
		/* Clear moving row/column in canvas*/
		clearCanvasRowColumn(dragParams.itemLabel,dragParams.type);
    }
    function dragLabelStop() { // also fired on click
		if (dragContext) { // do not set dragContext to null to prevent click action
			/* Stabilizing position */
			moveLabelAndCells(this,true);
			/* Delete temporary svg */
			for (let i=0; i<dragParams.svgCells.length; i++) {
				dragParams.svgCells[i].remove();
			}
			/* Redraw heat cells value lines */
			if (cellsValueTrace) {
				if (dragParams.movedLabelRanks.min > 1) {
					dragParams.movedLabelRanks.min--;
					moveCells(dragParams.labelList[dragParams.movedLabelRanks.min],dragParams.axis,0,true); // also redraw cells for lower-flanking label
				}
				if (dragParams.movedLabelRanks.max < dragParams.maxRank) {
					dragParams.movedLabelRanks.max++;
					moveCells(dragParams.labelList[dragParams.movedLabelRanks.max],dragParams.axis,0,true); // also redraw cells for upper-flanking label
				}
				for (let rank=dragParams.movedLabelRanks.min; rank<=dragParams.movedLabelRanks.max; rank++) {
					drawCellsValueTrace(dragParams.labelList[rank].heatCellList);
				}
			}
			/* Return new label-order */
			if (labelOnMove) {
				let labelList=(dragParams.type=='row')? rowList : columnList,
					labelIdList=[];
				for (let i=0; i<labelList.length; i++) {labelIdList.push(labelList[i].id);}
				labelOnMove(dragParams.type,labelIdList);
			}
		}
		dragParams={};
	}
    function dragLabelMove(dx,dy) {
		if (!dragContext) { // 1st move -> initialize dragParams
			dragLabelInit(this);
		}
		if (labelPopupText) {labelPopupText.remove(); labelPopupText=null;}
		dragContext=true;
		let dragLength;
		/* Moving label*/
		if (this.data('row')) { // row
			if (dy < 0) {dragLength=Math.max(labelPosRanges.minY-dragParams.startPosL,dy);}
			else {dragLength=Math.min(labelPosRanges.maxY-dragParams.startPosL,dy);}
			this.attr({y:dragParams.startPosL+dragLength});
		}
		else { // column
			if (dx < 0) {dragLength=Math.max(labelPosRanges.minX-dragParams.startPosL,dx);} // column text is rotated x->y & y->x
			else {dragLength=Math.min(labelPosRanges.maxX-dragParams.startPosL,dx);}
			// Move (+/-rotation)
			//this.rotate(-colLabelAngle,this.attr('x'),this.attr('y')); // unrotate
			let tr=this.transform().toString().replace(/r[^t]+/,''); // in case inside focusArea
			//this.transform(''); // unrotate
			this.attr(dragParams.axis,dragParams.startPosL+dragLength);
			//this.rotate(colLabelAngle,this.attr('x'),this.attr('y')); // rerotate
			this.transform(tr+'r'+colLabelAngle+','+this.attr('x')+','+this.attr('y')); // rerotate
		}
		/* Moving svg cells */
		let dragDelta=dragLength-dragParams.dragPrevLength;
		for (let i=0; i<dragParams.svgCells.length; i++) {
			dragParams.svgCells[i].attr(dragParams.axis,dragParams.svgCells[i].attr(dragParams.axis)+dragDelta);
		}
		/* Moving cell selection (if any) */
		if (dragParams.svgCells.length > 1) { // tempory svg cells
			let shiftX=0,shiftY=0;
			if (dragParams.axis=='x') {shiftX=this.attr(dragParams.axis)-dragParams.emptyPosL;}
			else {shiftY=this.attr(dragParams.axis)-dragParams.emptyPosL;}
			for (let i=0; i<this.data(dragParams.type).heatCellList.length; i++) {
				let hCell=this.data(dragParams.type).heatCellList[i];
				if (hCell.isSelectedData) {
					hCell.isSelectedData.remove();
					drawCellSelection(hCell,shiftX,shiftY);
				}
			}
		}
		/* moving annotation cells */
		//var shift=dragParams.emptyPosL-this.attr(dragParams.axis);
		//for (let setName in dragParams.itemLabel.annotations) {
        Object.keys(dragParams.itemLabel.annotations).forEach(function(setName) {
		    let aC=dragParams.itemLabel.annotations[setName].aspect;
			aC.attr(dragParams.axis,aC.attr(dragParams.axis)+dragDelta);
		});
		/* Checking if dragging over another label */
		let movedLabelRank=null;
		if (dragDelta < 0) { // dragging up/left
			if (dragParams.itemLabel.rank > 1) {
				let overLabel=dragParams.labelList[dragParams.itemLabel.rank-1];
				if (this.attr(dragParams.axis)-overLabel.aspect.attr(dragParams.axis) < maxOverlap[dragParams.type]) {
					clearCanvasRowColumn(overLabel,dragParams.type);
					movedLabelRank=overLabel.rank; // record before moveLabelAndCells()
					moveLabelAndCells(overLabel.aspect); // false => move canvas cells
				}
			}
		}
		else { // dragging down/right
			if (dragParams.itemLabel.rank < dragParams.maxRank) { // max rank
				let overLabel=dragParams.labelList[dragParams.itemLabel.rank+1];
				if (overLabel.aspect.attr(dragParams.axis)-this.attr(dragParams.axis) < maxOverlap[dragParams.type]) { //-2 for safety
					clearCanvasRowColumn(overLabel,dragParams.type);
					movedLabelRank=overLabel.rank; // record before moveLabelAndCells()
					moveLabelAndCells(overLabel.aspect); // false => move canvas cells
				}
			}
		}
		if (movedLabelRank) {
			dragParams.movedLabelRanks.min=Math.min(dragParams.movedLabelRanks.min,movedLabelRank);
			dragParams.movedLabelRanks.max=Math.max(dragParams.movedLabelRanks.max,movedLabelRank);
		}
		dragParams.dragPrevLength=dragLength;
    }
	function clearCanvasRowColumn(label,type) {
		let setX,setY,setW,setH; // row/column of cells selected for drag
		if (type=='row') {
			setX=0; //cellX0;
			setY=hCellH*(label.rank-1);
			setW=hCellW*label.heatCellList.length;
			setH=hCellH;
		}
		else { // column
			setX=hCellW*(label.rank-1); //cellX0+
			setY=0; //cellY0;
			setW=hCellW;
			setH=hCellH*label.heatCellList.length;
		}
		mapContext.fillStyle='#FFFFFF';
		mapContext.fillRect(setX,setY,setW,setH);
	}
    function moveLabelAndCells(movedLabel,endDrag) {
		let prevPosL=movedLabel.attr(dragParams.axis);
		/* moving canvas cells */
		let shift=dragParams.emptyPosL-prevPosL;
		moveCells(movedLabel.data(dragParams.type),dragParams.axis,shift,endDrag);
		/* moving label */
		let dragList,otherList,dragTypeIdx;
		if (dragParams.type=='row') {
			movedLabel.attr({y:dragParams.emptyPosL});
			dragList=rowList;
			otherList=columnList;
			dragTypeIdx='rowIdx';
		}
		else { // column
			//movedLabel.rotate(-colLabelAngle,movedLabel.attr('x'),movedLabel.attr('y')); // unrotate
			movedLabel.attr({x:dragParams.emptyPosL});
			//movedLabel.rotate(colLabelAngle,movedLabel.attr('x'),movedLabel.attr('y')); // rerotate
			dragList=columnList;
			otherList=rowList;
			dragTypeIdx='columnIdx';
		}
		/* switching ranks & indexes */
		let overRank=movedLabel.data(dragParams.type).rank,
			itemRank=dragParams.itemLabel.rank;
		if (overRank != itemRank) { // not end of dragging
			movedLabel.data(dragParams.type).rank=itemRank;
			dragParams.itemLabel.rank=overRank;
			dragParams.labelList[itemRank]=movedLabel.data(dragParams.type);
			dragParams.labelList[overRank]=dragParams.itemLabel;
			//indexes
			let overCellList=dragList[overRank-1];
			dragList[overRank-1]=dragList[itemRank-1];
			dragList[itemRank-1]=overCellList;
			for (let i=0; i<otherList.length; i++) {
				let overCell=otherList[i].heatCellList[overRank-1];
				otherList[i].heatCellList[overRank-1]=otherList[i].heatCellList[itemRank-1];
				otherList[i].heatCellList[itemRank-1]=overCell;
				otherList[i].heatCellList[overRank-1][dragTypeIdx]=overRank-1;
				otherList[i].heatCellList[itemRank-1][dragTypeIdx]=itemRank-1;
			}
		}
		/* Label rotation & translation if inside focusArea */
		let overIdx=itemRank-1; // label were switched or end of drag
		if (dragParams.type=='row') {
			let trF='';
			if (focusArea.data('status')=='on' && overIdx >= focusArea.data('labelIdx')[0] && overIdx <= focusArea.data('labelIdx')[1]) { // inside focus zone
				let focusShift=(rowLabelLocation=='left')? focusArea.attr('x')-cellX0 : focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
				trF='t'+focusShift+',0';
				movedLabel.toFront(); // shift & move above bgRow
			}
			movedLabel.transform(trF); // also when empty in case moved outside focus
		}
		else {
			let trF='';
			if (focusArea.data('status')=='on' && overIdx >= focusArea.data('labelIdx')[2] && overIdx <= focusArea.data('labelIdx')[3]) { // inside focus zone
				let focusShift=(colLabelLocation=='top')? focusArea.attr('y')-cellY0 : focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
				trF='t0,'+focusShift;
				movedLabel.toFront(); // move above bgColumn
			}
			movedLabel.transform(trF+'r'+colLabelAngle+','+movedLabel.attr('x')+','+movedLabel.attr('y'));
		}

		/* storing new empty pos */
		dragParams.emptyPosL=prevPosL;
    }

    function moveCells(label,axis,shift,endDrag) { // axis: x (column) or y (row)!
		let cList=label.heatCellList,
			shiftX,shiftY;
		if (axis=='y') {
			//label=rowList[cList[0].rowIdx];
			shiftX=0;
			shiftY=(endDrag)? 0 : shift;
		}
		else { // x
			//label=columnList[cList[0].columnIdx];
			shiftX=(endDrag)? 0 : shift;
			shiftY=0;
		}
		/* moving heat cells */
		for (let i=0; i < cList.length; i++) {
			//canvas
			mapContext.fillStyle=cList[i].color;
			mapContext.fillRect((hCellW*cList[i].columnIdx)+shiftX,(hCellH*cList[i].rowIdx)+shiftY,hCellW,hCellH);
			if (cList[i].isFlagged && cellFlagVisible) { // cell flag
				drawCellFlag((hCellW*cList[i].columnIdx)+shiftX,(hCellH*cList[i].rowIdx)+shiftY);
			}
			//svg
			if (cList[i].isSelectedData) {
				cList[i].isSelectedData.remove();
				drawCellSelection(cList[i],shiftX,shiftY);
			}
		}

		/* moving annotation cells */
		//for (let setName in label.annotations) {
        Object.keys(label.annotations).forEach(function(setName) {
		    let aC=label.annotations[setName].aspect;
			aC.attr(axis,aC.attr(axis)+shift);
		});
    }

	this.renameLabel=function(type,labelID,newName) {
		let labelList=(type=='row')? rowList : columnList;
		for (let i=0; i<labelList.length; i++) {
			if (labelList[i].id==labelID) {
				labelList[i].label=newName;
				labelList[i].aspect.attr('text',newName);
				break;
			}
		}
	};


	/**********************************************************/
    /************** Drag-drawing zoom area events *************/
    /**********************************************************/
	function dragFocusStart(x,y,e) { // x,y are absolute to document // *public*
		if (focusArea.data('status') == 'on') return; // because triggered by simple click (prevents reset of focusArea on simple click)
		//var mouseInChart=getMousePositionInChart(e);
//debug('ctrl='+idvis.lib.modKeyPressed+' : '+x+'('+mousePos[0]+'), '+y+'('+mousePos[1]+')');
		let mousePos=idvis.getMousePositionInElement(e),
			mouseX=mousePos[0]+cellX0,mouseY=mousePos[1]+cellY0;
		/*
		if (e.offsetX) { // IE
			if (e.layerX) { // IE 9+
				mouseX=e.offsetX;
				mouseY=e.offsetY;
			}
			else { // IE 5-8
				mouseX=e.offsetX+cellX0;
				mouseY=e.offsetY+cellY0;
			}
		}
		else {
			mouseX=e.layerX;
			mouseY=e.layerY;
		}
		*/
//		dragContext=true; // ?
		hCellsSensor.attr({cursor:'crosshair'});

		//if (focusArea.data('status') == 'on') {
		//	focusArea.attr({width:0,height:0});
		//}
		focusArea.data('close')[0].attr({fill:cellHighlightColor});
		if (focusArea.data('selRows')) {focusArea.data('selRows')[0].attr({fill:cellHighlightColor});}
		if (focusArea.data('selColumns')) {focusArea.data('selColumns')[0].attr({fill:cellHighlightColor});}
		focusArea.data({startX:mouseX,startY:mouseY}) //,status:'extend'
		.attr({x:mouseX,y:mouseY,width:0,height:0,stroke:cellHighlightColor,fill:cellHighlightColor})
		.toFront()
		.show();
//console.log('DRAG_0');
	}
	function dragFocusMove(dx,dy,x,y,e) { // *public*
		if (focusArea.data('status') == 'on') {return;} // prevents extension of already drawn focusArea
        dragContext=true; // allow to distinguish a drag from a click!
		if (dx > 0) {
			focusArea.attr({width:dx});
			if (focusArea.attr('x')+focusArea.attr('width') > mapW+cellX0) {focusArea.attr({width:mapW-focusArea.attr('x')+cellX0});}
		}
		else {
			focusArea.attr({x:focusArea.data('startX')+dx,width:-dx});
			if (focusArea.attr('x') < cellX0) {
				let extra=cellX0 - focusArea.attr('x');
				focusArea.attr({x:cellX0,width:focusArea.attr('width')-extra});
			}
		}
		if (dy > 0) {
			focusArea.attr({height:dy});
			if (focusArea.attr('y')+focusArea.attr('height') > mapH+cellY0) {focusArea.attr({height:mapH-focusArea.attr('y')+cellY0});}
		}
		else {
			focusArea.attr({y:focusArea.data('startY')+dy,height:-dy});
			if (focusArea.attr('y') < cellY0) {
				let extra=cellY0 - focusArea.attr('y');
				focusArea.attr({y:cellY0,height:focusArea.attr('height')-extra});
			}
		}
	}
	function dragFocusEnd() { // *public*
        
        /* Click not drag ? */      
        if (!dragContext) { // click NOT drag event
//console.log('CLICK2',idvis.modKeyPressed);
            if (cellOnClick) {
                if (idvis.modKeyPressed && !noCellExclusion) { excludeCell(lastHighlightedCell); } else if (lastHighlightedCell) { selectCell(lastHighlightedCell); } //else {console.log('OTHER');}
            }
            else if (!noCellExclusion && idvis.modKeyPressed) { excludeCell(lastHighlightedCell); }
        }

		if (focusArea.data('status') == 'on') {return;} // because triggered by simple click
		dragContext=false;
		if (focusArea.attr('width') >= hCellW && focusArea.attr('height') >= hCellH) { // big enough => end extension
			focusArea.data({status:'on'});
			/* Find matching row & column label ranges */
			let colIdx1=Math.max(0,Math.round(((focusArea.attr('x')-cellX0)/hCellW)-0.5)),
				colIdx2=Math.min(columnList.length-1,Math.round(((focusArea.attr('x')+focusArea.attr('width')-cellX0)/hCellW)-0.5)),
				rowIdx1=Math.max(0,Math.round(((focusArea.attr('y')-cellY0)/hCellH)-0.5)),
				rowIdx2=Math.min(rowList.length-1,Math.round(((focusArea.attr('y')+focusArea.attr('height')-cellY0)/hCellH)-0.5));
//console.log(colIdx1,colIdx2,rowIdx1,rowIdx2);
			/* Adjust drag area */
			focusArea.attr({x:cellX0+hCellW*colIdx1,
						  y:cellY0+1+hCellH*rowIdx1,
						  width:hCellW*(colIdx2-colIdx1+1),
						  height:hCellH*(rowIdx2-rowIdx1+1)
						})
			.data({labelIdx:[rowIdx1,rowIdx2,colIdx1,colIdx2]});

			/* Reset cell sensor */
			hCellsSensor.toFront();
			/* Display focusArea control(s) */
			//let numRowAnnot=0,numColAnnot=0;
			//for (let sName in annotationList.row) {numRowAnnot++;}
			//for (let sName in annotationList.column) {numColAnnot++;}
            let numRowAnnot=Object.keys(annotationList.row).length,
                numColAnnot=Object.keys(annotationList.column).length,
			    annotSpaceX=numRowAnnot*annotationCellSize.row,
				annotSpaceY=numColAnnot*annotationCellSize.column,
				cX=(rowLabelLocation=='left')? focusArea.attr('x')+focusArea.attr('width')+annotSpaceX-6 : focusArea.attr('x')-annotSpaceX-28,
				setY=focusArea.attr('y')-6;
			focusArea.data('close').transform('T'+cX+','+setY).show().toFront();
			if (focusArea.data('selRows')) {focusArea.data('selRows').transform('T'+cX+','+setY).show().toFront();}
			if (focusArea.data('selColumns')) {focusArea.data('selColumns').transform('T'+cX+','+setY).show().toFront();}

			/* Move labels & annotations around focus area */
			//Rows
			let bgW=0;
			for (let i=rowIdx1; i<=rowIdx2; i++) {bgW=Math.max(bgW,rowList[i].aspect.getBBox().width);}
			bgW+=10;
			if (rowLabelLocation=='left') {
				if (hCellH >= 5) focusArea.data('bgRow').attr({x:focusArea.attr('x')-bgW,y:focusArea.attr('y')-2,width:bgW,height:focusArea.attr('height')+4});
				focusArea.data('bgRow2').attr({x:focusArea.attr('x')+focusArea.attr('width')+2,y:focusArea.attr('y'),width:annotSpaceX,height:focusArea.attr('height')});
			}
			else {
				focusArea.data('bgRow').attr({x:focusArea.attr('x')+focusArea.attr('width'),y:focusArea.attr('y')-4,width:bgW,height:focusArea.attr('height')+8});
				focusArea.data('bgRow2').attr({x:focusArea.attr('x')-annotSpaceX-2,y:focusArea.attr('y'),width:annotSpaceX,height:focusArea.attr('height')});
			}
			if (hCellH >= 5) focusArea.data('bgRow').show().toFront();
			focusArea.data('bgRow2').show().toFront();
			let shiftX,shiftX2;
			if (rowLabelLocation=='left') {
				shiftX=focusArea.attr('x')-cellX0;
				shiftX2=focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
			}
			else {
				shiftX=focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
				shiftX2=focusArea.attr('x')-cellX0;
			}
			for (let i=rowIdx1; i<=rowIdx2; i++) {
				rowList[i].aspect.toFront().transform('t'+shiftX+',0');
				/* Moving annotation cells */
				for (let setName in rowList[i].annotations) {
					rowList[i].annotations[setName].aspect.toFront().transform('t'+shiftX2+',0');
				}
			}
			//Columns
			let bgH=0;
			for (let i=colIdx1; i<=colIdx2; i++) {bgH=Math.max(bgH,columnList[i].aspect.getBBox().height);}
			bgH+=10;
			if (colLabelLocation=='top') {
				let shiftD=(colLabelOrientation=='diagonal')? columnList[colIdx2].aspect.getBBox().width : 0;
				focusArea.data('bgColumn').attr({x:focusArea.attr('x')-2,y:focusArea.attr('y')-bgH,width:focusArea.attr('width')+shiftD+4,height:bgH});
				focusArea.data('bgColumn2').attr({x:focusArea.attr('x'),y:focusArea.attr('y')+focusArea.attr('height')+2,width:focusArea.attr('width'),height:annotSpaceY});
			}
			else {
				let shiftD=(colLabelOrientation=='diagonal')? columnList[colIdx1].aspect.getBBox().width : 0;
				focusArea.data('bgColumn').attr({x:focusArea.attr('x')-shiftD-2,y:focusArea.attr('y')+focusArea.attr('height'),width:focusArea.attr('width')+shiftD+4,height:bgH});
				focusArea.data('bgColumn2').attr({x:focusArea.attr('x'),y:focusArea.attr('y')-annotSpaceY-2,width:focusArea.attr('width'),height:annotSpaceY});
			}
			focusArea.data('bgColumn').show().toFront();
			focusArea.data('bgColumn2').show().toFront();
			let shiftY,shiftY2;
			if (colLabelLocation=='top') {
				shiftY=focusArea.attr('y')-cellY0;
				shiftY2=focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
			}
			else {
				shiftY=focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
				shiftY2=focusArea.attr('y')-cellY0;
			}
			for (let i=colIdx1; i<=colIdx2; i++) {
				let tr0=columnList[i].aspect.transform();
				columnList[i].aspect.toFront().transform('t0,'+shiftY+tr0);
//console.log('DRAG: '+columnList[i].aspect.transform());  //.toString()
				/* moving annotation cells */
				for (let setName in columnList[i].annotations) {
					columnList[i].annotations[setName].aspect.toFront().transform('t0,'+shiftY2);
				}
			}
		}
		else { // too small => hide
			focusArea.hide()
			.data({status:'off'})
			.attr({width:0,height:0});
		}
		hCellsSensor.attr({cursor:'auto'});
	}
	function closeFocusArea() {
		if (!hasCells || focusArea.data('status')=='off') return; // 1D or already closed
		focusArea.hide().data({status:'off'}).data('close').hide();
		if (focusArea.data('selRows')) {focusArea.data('selRows').hide();}
		if (focusArea.data('selColumns')) {focusArea.data('selColumns').hide();}
		focusArea.data('bgRow').hide();
		focusArea.data('bgRow2').hide();
		focusArea.data('bgColumn').hide();
		focusArea.data('bgColumn2').hide();
		// Rows: delete transformation
		for (let i=focusArea.data('labelIdx')[0]; i<=focusArea.data('labelIdx')[1]; i++) {
			rowList[i].aspect.transform('');
			for (let setName in rowList[i].annotations) { // annotation cells
				rowList[i].annotations[setName].aspect.transform('');
			}
		}
		//Columns: remove translation, keep original rotation
		for (let i=focusArea.data('labelIdx')[2]; i<=focusArea.data('labelIdx')[3]; i++) {
			//var tr1=columnList[i].aspect.transform().toString().replace(/t[^r]+/,'');
			//tr1=tr1.replace(/,.+/,','+columnList[i].aspect.attr('x')+','+columnList[i].aspect.attr('y')); // transform() has forgotten x,y coordinate since >raphael 2.12!!!!!
//console.log('END: '+columnList[i].aspect.transform(),'  ->  ',tr1,' -> ',tr1);  //.toString()
			columnList[i].aspect.transform('r'+colLabelAngle+','+columnList[i].aspect.attr('x')+','+columnList[i].aspect.attr('y')); // rerotate text
//console.log(''+columnList[i].aspect.transform());
			for (var setName in columnList[i].annotations) { // annotation cells
				columnList[i].annotations[setName].aspect.transform('');
			}
		}
	}

	/**********************************************************/
    /******************** Annotation events *******************/
    /**********************************************************/
    /*** Adding Annotations ***/
    this.addAnnotation=function(type,setName,setData,matchPattern) { // (matchPattern: optional. use ### in matchPattern as point id)
		if (type=='row' && rowList.length===0) {return true;} // no rows => 1D-vertical clustering BUT 'true'=silent validation
		//Check if exists already
		if (annotationList[type][setName]) {
			alert('Annotation \''+setName+'\' already displayed!');
			return false;
		}
		let startShift=annotationCellSize[type]+2,
			minSpace=(branchList[type])? 50 : -50;
		if (treeSpace[type]-startShift < minSpace) {
			alert('Sorry: No space left. Please delete 1 annotation first.');
			return false;
		}
		closeFocusArea(); // to prevent dealing with its update
		treeSpace[type]-=startShift;
		let labelList,
			annotStart=treeStart[type];
		if (type=='row') {
			labelList=rowList;
			annotStart=(rowLabelLocation=='left')? treeStart[type] : treeStart[type]-startShift;
			treeStart.row+=(rowLabelLocation=='left')? startShift : -startShift;
		}
		else {
			labelList=columnList;
			annotStart=(colLabelLocation=='top')? treeStart[type] : treeStart[type]-startShift;
			treeStart.column+=(colLabelLocation=='top')? startShift : -startShift;
		}
		clearTree(type);
		drawTree(type);

	    let labelIdList={}; // id -> Label (list of available ids)
	    for (let i=0; i<labelList.length; i++) {
			labelIdList[labelList[i].id]=labelList[i];
	    }
		let annotationSet=new AnnotationSet(type,setName),
			l; // set label SVG object
		if (type=='row') { // Row annot name written alongside the column labels
			let lX=annotStart+Math.round(annotationCellSize[type]/2),
				lY,lAnchor;
			//if (colLabelLocation=='top') {
				lY=cellY0-5; // always top not to interfere with legends (PP 02/10/16)
				lAnchor='start';
			//}
			//else {
			//	lY=cellY0+mapH+5;
			//	lAnchor='end';
			//}
			//var dispName=(setName.length <= maxColLabelSize)? setName : idvis.lib.shortenText(setName,maxColLabelSize);
			//l=HM.paper.text(lX,lY,dispName).attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':lAnchor}).transform('r'+colLabelAngle+','+lX+','+lY);
			l=HM.paper.text(lX,lY,setName).attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':lAnchor});
			let fontSize=colLabelFontSize;
			while (l.getBBox().width*colLabelRatio > colLabelSpace) {
				if (fontSize <= 6) {
					let ratio=colLabelSpace/l.getBBox().width*colLabelRatio,
						dispName=idvis.lib.shortenText(setName,Math.round(setName.length*ratio));
					l.attr({'text':dispName});
					break;
				}
				fontSize-=2;
				l.attr({'font-size':fontSize});
			}
			let angle=(rowLabelLocation=='left' && colLabelLocation=='top')? colLabelAngle : -90; // -90 no diagonal to prevent colision with tree/annot cell
			l.transform('r'+angle+','+lX+','+lY);
		}
		else { // Column annot name written alongside the row labels
			let lY=annotStart+Math.round(annotationCellSize[type]/2),
				lX,lAnchor;
			if (rowLabelLocation=='left') {
				lX=cellX0-5;
				lAnchor='end';
			}
			else {
				lX=cellX0+mapW+5;
				lAnchor='start';
			}
			//var dispName=(setName.length <= maxRowLabelSize)? setName : idvis.lib.shortenText(setName,maxRowLabelSize);
			//l=HM.paper.text(lX,lY,dispName).attr({'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':lAnchor});
			l=HM.paper.text(lX,lY,setName).attr({'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':lAnchor});
			let fontSize=rowLabelFontSize;
			while (l.getBBox().width > rowLabelSpace) {
				if (fontSize <= 6) {
					let ratio=rowLabelSpace/l.getBBox().width,
						dispName=idvis.lib.shortenText(setName,Math.round(setName.length*ratio));
					l.attr({'text':dispName});
					break;
				}
				fontSize-=2;
				l.attr({'font-size':fontSize});
			}
		}
		l.data('set',annotationSet)
		//.hover(function(){this.attr({fill:'#F00',cursor:'pointer'});},function(){this.attr({fill:'#000',cursor:'auto'});})
		.hover(function(){setLabelEmphasis(this,'on');},function(){setLabelEmphasis(this,'off');})
		.dblclick(function(){if (confirm('Delete annotation \''+setName+'\'?')) {deleteAnnotation(this);}});

		annotationSet.aspect=l;
		let usedTargetIds={},
			colorIdx=0; // 0++ -> 1: 0 not used for annotations. Index reset for each annot (PP 29/0516)
		
        //Color:Scan for matching previously used annot name => use same color   
        var reservedColor={};
        for (let s=0; s<setData.length; s++) {
            let [annot,targetData,annotColor]=setData[s]; // name,targetStrg,color(optional)
            if (!annotColor && annotNameToColor[type][annot]) { // already seen
                reservedColor[annotNameToColor[type][annot]]=1;
            }
        }

        for (let s=0; s<setData.length; s++) {
			let [annot,targetData,annotColor]=setData[s]; // name,targetStrg,color(optional)
			//Color
			if (!annotColor) { // Choose color from default list
                if (annotNameToColor[type][annot]) { // already seen
                    annotColor=annotNameToColor[type][annot];
                }
                else {
                    var firstColorIdx;
                    var firstLoop=true;
                    while (true) { // loop through all colors to find a free one. If none: use next in list anyway
                        colorIdx++;    
                        if (colorIdx >= colorList.length) colorIdx=1; // starts at 1 (0: black!)
                        if (firstLoop) {
                            firstColorIdx=colorIdx;
                            firstLoop=false;
                        }
                        else if (colorIdx == firstColorIdx) {break;}
                        if (!reservedColor[colorList[colorIdx]]) break;
                    }
                    annotColor=colorList[colorIdx];
                    annotNameToColor[type][annot]=annotColor;
                }
            }
            //Data
			let targetLabels=targetData.split(',');
			for (let i=0; i<targetLabels.length; i++) {
				//if (usedTargetIds[targetLabels[i]]) continue; // ID/label already used by this annotation
				if (matchPattern) {
					/* loose ID match => loop through all label Ids */
					let matchRegExp=new RegExp(matchPattern.replace('###',targetLabels[i]));
					for (let labelID in labelIdList) {
						//if (usedTargetIds[labelID]) continue; // ID/label already used by this annotation
						if (matchRegExp.exec(labelID)) {
							if (usedTargetIds[labelID]) { // ID/label already used by this annotation
								let annotCell=usedTargetIds[labelID];
								annotCell.name+='\n'+annot;
								annotCell.color='#000000';
								annotCell.aspect.attr({fill:'#000000'});
							}
							else {
								let label=labelIdList[labelID];
								usedTargetIds[labelID]=addAnnotationCell(type,label,annotationSet,setName,annot,annotStart,annotationCellSize[type],annotColor);
								//usedTargetIds[labelID]=1;
							}
						}
					}
				}
				else {
					/* ID equality (single match) */
					//if (usedTargetIds[targetLabels[i]]) continue; // ID/label already used by this annotation
					if (usedTargetIds[targetLabels[i]]) { // ID/label already used by this annotation
						let annotCell=usedTargetIds[targetLabels[i]];
						annotCell.name+='\n'+annot;
						annotCell.color='#000000';
						annotCell.aspect.attr({fill:'#000000'});
					}
					else if (labelIdList[targetLabels[i]]) { // make sure ID is in known
						usedTargetIds[targetLabels[i]]=addAnnotationCell(type,labelIdList[targetLabels[i]],annotationSet,setName,annot,annotStart,annotationCellSize[type],annotColor);
						//usedTargetIds[targetLabels[i]]=1;
					}
				}
				//
				//var annotCell=new AnnotationCell(annotationSet,annot,label,annotColor);
				///* Drawing annotation cell */
				//var x,y,w,h;
				//if (type=='row') {
				//	x=(rowLabelLocation=='left')? annotStart : annotStart+2;
				//	y=cellY0+hCellH*(label.rank-1);
				//	w=cellSize;
				//	h=hCellH;
				//}
				//else {
				//	x=cellX0+hCellW*(label.rank-1);
				//	y=(colLabelLocation=='top')? annotStart : annotStart+2;
				//	w=hCellW;
				//	h=cellSize;
				//}
				//var aC=HM.paper.rect(x,y,w,h,3).attr({fill:annotColor,stroke:'none','fill-opacity':0.5}) //stroke:colorList[colorIdx],'stroke-width':0,'stroke-opacity':0.2
				//.data('annotCell',annotCell)
				//.hover(function(e){highlightAnnotCell(this,'on',e)},function(){highlightAnnotCell(this,'off')});
				//annotCell.aspect=aC;
				//label.annotations[setName]=annotCell;
				//annotationSet.cells.push(annotCell);
				//usedTargetIds[targetLabels[i]]=1;
			}
		}
		let maxRank=0;
		for (let sName in annotationList[type]) {maxRank=Math.max(maxRank,annotationList[type][sName].rank);}
		annotationSet.rank=maxRank+1;
		annotationList[type][setName]=annotationSet;
		updateAnnotationLegends();
		return true;
	};

	function addAnnotationCell(type,label,annotationSet,setName,annot,annotStart,cellSize,annotColor) {
		let annotCell=new AnnotationCell(annotationSet,annot,label,annotColor);
		/* Drawing annotation cell */
		let x,y,w,h;
		if (type=='row') {
			x=(rowLabelLocation=='left')? annotStart : annotStart+2;
			y=cellY0+hCellH*(label.rank-1);
			w=cellSize;
			h=hCellH;
		}
		else {
			x=cellX0+hCellW*(label.rank-1);
			y=(colLabelLocation=='top')? annotStart : annotStart+2;
			w=hCellW;
			h=cellSize;
		}
		let aC=HM.paper.rect(x,y,w,h,3).attr({fill:annotColor,stroke:'none'}) //,'fill-opacity':0.75 stroke:colorList[colorIdx],'stroke-width':0,'stroke-opacity':0.2
		.data('annotCell',annotCell)
		.hover(function(e){highlightAnnotCell(this,'on',e);},function(){highlightAnnotCell(this,'off');});
		annotCell.aspect=aC;
		label.annotations[setName]=annotCell;
		annotationSet.cells.push(annotCell);
		return annotCell;
	}
	function deleteAnnotation(l) {
		let annotationSet=l.data('set'),
			setName=annotationSet.name,
			type=annotationSet.type;
		if (annotationOnDelete) {
			let res=annotationOnDelete(type,setName);
			if (!res) {
				alert('An external error has occured. Cannot delete \''+setName+'\'');
				return;
			}
		}
		let annotRank=annotationSet.rank;
		// Delete annotation cells
		for (let i=0; i<annotationSet.cells.length; i++) {
			annotCell=annotationSet.cells[i];
			annotCell.aspect.removeData('annotCell');
			annotCell.aspect.remove();
			annotCell.label.annotations[setName]=null;
			delete annotCell.label.annotations[setName];
			annotCell=null;
		}
		annotationSet.cells=null;
		l.removeData('set');
		l.remove();
		annotationSet.aspect=null;
		annotationSet=null;
		annotationList[type][setName]=null;
		delete annotationList[type][setName];
		//Shift down all following annotation sets
		let cellSize=(type=='row')? Math.min(hCellW,20) : (hasCells)? Math.min(hCellH,20) : 20, // 1/2D for type=column (same in AddAnnotation)
			dx=0,dy=0;
		if (type=='row') {dx=(rowLabelLocation=='left')? -cellSize-2 : cellSize+2;}
		else {dy=(colLabelLocation=='top')? -cellSize-2 : cellSize+2;}
		for (let sName in annotationList[type]) {
			if (annotationList[type][sName].rank > annotRank) {
				annotationList[type][sName].rank--;
				let l2=annotationList[type][sName].aspect;
				if (type=='row') {
					//l2.rotate(-colLabelAngle,l2.attr('x'),l2.attr('y')); // unrotate
					l2.transform(''); // unrotate
					l2.attr({x:l2.attr('x')+dx});
					//l2.rotate(colLabelAngle,l2.attr('x'),l2.attr('y')); // rerotate
					l2.transform('r'+colLabelAngle+','+l2.attr('x')+','+l2.attr('y')); // rerotate
				}
				else {l2.attr({y:l2.attr('y')+dy});}
				for (let i=0; i<annotationList[type][sName].cells.length; i++) {
					let aC=annotationList[type][sName].cells[i].aspect;
					aC.attr({x:aC.attr('x')+dx,y:aC.attr('y')+dy});
				}
			}
		}
		treeStart.row+=dx;
		treeSpace.row+=Math.abs(dx);
		treeStart.column+=dy;
		treeSpace.column+=Math.abs(dy);
		clearTree(type);
		drawTree(type);
		closeFocusArea(); // to prevent dealing with its update
		updateAnnotationLegends();
	}

	function updateAnnotationLegends() {
		annotationLegends.remove(); // clear previous legend if any
		let posX=(branchList.row || rowLabelLocation=='right')? paperWidth0 : treeStart.row + 30, //cellX0+mapW+20;
			posY=(hasCells)? cellY0 : 30, //(hasCells)? cellY0+mapH+20 : 20;
			types=(hasCells)? ['column','row'] : ['column'];
		for (let t=0; t<2; t++) {
			var type=types[t];
			var numAnnot=0;
			posY+=10;
			for (let setName in annotationList[type]) {
				numAnnot++;
				if (numAnnot==1 && hasCells) {
					annotationLegends.push(HM.paper.text(posX,posY,HM.entities[type]).attr({'font-size':14,'font-weight':'bold','text-anchor':'start'}));
					posY+=15;
				}
				annotationLegends.push(HM.paper.text(posX+10,posY,setName+':').attr({'font-size':12,'font-weight':'bold','text-anchor':'start'}));
				posY+=15;
				let usedCells={};
				for (let i=0; i<annotationList[type][setName].cells.length; i++) {
					let annotName=annotationList[type][setName].cells[i].name;
					if (annotName.match(/\n/)) {annotName='*Multi-label*';}
					else {
						let annotNameInfo=annotName.split('##');
						annotName=annotNameInfo[0];
					}
					if (usedCells[annotName]) {continue;}
					let annotColor=annotationList[type][setName].cells[i].color;
					annotationLegends.push(HM.paper.rect(posX+15,posY-5,10,10).attr({fill:annotColor,stroke:annotColor}));
					annotationLegends.push(HM.paper.text(posX+30,posY,annotName).attr({'font-size':12,'text-anchor':'start'}));
					usedCells[annotName]=1;
					posY+=15;
				}
			}
		}
//console.log(annotationLegends.getBBox());
		/*Readjust paper size if necessary */
		let legend=annotationLegends.getBBox(),
            prevPaperWidth=paperWidth,
            prevPaperHeight=paperHeight;
		paperWidth=Math.max(paperWidth0,legend.x2+20);
		paperHeight=Math.max(paperHeight0,legend.y2+20);
		HM.paper.setSize(paperWidth,paperHeight);
		let pathArray=backgroundPath.attr('path');
		pathArray[1][1]=paperWidth;
		pathArray[2][1]=paperWidth; pathArray[2][2]=paperHeight;
		pathArray[3][2]=paperHeight;
		backgroundPath.attr({path:pathArray});
		background.attr({width:paperWidth,height:paperHeight});
        let dWidth=paperWidth-prevPaperWidth,
            dHeight=paperHeight-prevPaperHeight;
        mainDiv.style.width=parseInt(mainDiv.style.width) + dWidth + 'px'; // extend main outer div width
        mainDiv.style.height=parseInt(mainDiv.style.height) + dHeight + 'px'; // extend main outer div height
	}

	function highlightAnnotCell(aC,action,e) {
		let annotCell=aC.data('annotCell');
		if (action=='on') {
			aC.attr({'stroke-width':2,stroke:cellHighlightColor});
			let x=aC.attr('x')+Math.round(aC.attr('width')/2),
				y=aC.attr('y'),
//console.log(aC.transform());
				textStrg='';
			if (e.ctrlKey || e.shiftKey) {
				textStrg+=HM.entities[annotCell.label.type];
				textStrg+=': '+annotCell.label.label+'\nSet: '+annotCell.set.name+'\n';
			}
			let multiAnnotInfo=annotCell.name.split('\n'); // in case multi annot on same cell
			for (let i=0; i<multiAnnotInfo.length; i++) {
				let annotNameInfo=multiAnnotInfo[i].split('##'); // e.g. 'annotX##1' -> annotX + 1
				if (i>0) {textStrg+='\n';}
				textStrg+=annotNameInfo[0];
			}
			drawCellPopup(x,y,textStrg,'#000','#FFF','#000');
			cellPopupText.transform(aC.transform()); // in case inside focusArea
		}
		else {
			aC.attr({'stroke-width':0,stroke:'none'});
			cellPopupText.remove();
		}
	}

	/**********************************************************/
    /********************** Tree events ***********************/
    /**********************************************************/
    function clearTree(type) {
		if (!branchList[type]) return; // no tree
		for (let i=0; i<branchList[type].length; i++) {
			let branch=branchList[type][i];
			branch.aspect[0].remove();
			branch.aspect[1].unclick(); // circle
			branch.aspect[1].removeData();
			branch.aspect[1].remove();
			branch.aspect.length=0;
		}
	}
	function clearTreePopup(type,tStrg) {
		if (treePopup[type][0].data('trString')==tStrg) { // Popup was not reset to new position by tree node click (.attr('transform')=[["t",x,y]] NOT string AND does not work in IE8!)
			treePopup[type].animate({'fill-opacity':0,'stroke-opacity':0},300,'linear',
									function() {
										treePopup[type].hide().transform('');
										treePopup[type][0].attr({'fill-opacity':0.7,'stroke-opacity':0.7}); // box
										treePopup[type][1].attr({'fill-opacity':1}); // text
									}
									);
		}
	}

	function drawTree(type) { // row or column
		if (!branchList[type]) return; // no tree
		let leaves={leaf1:{},leaf2:{}};
		for (let i=0; i<branchList[type].length; i++) {
//console.log(i);
			let branch=branchList[type][i];
			if (type=='row') { // row tree
				for (let leaf in leaves) {
					if (branch[leaf].heatCellList) { // label
						leaves[leaf].x=treeStart[type];
						leaves[leaf].y=branch[leaf].aspect.attr('y');
					}
					else { // branch
						branch[leaf].parent=branch;
						leaves[leaf].x=branch[leaf].connectionPos;
						leaves[leaf].y=branch[leaf].centerPos;
					}
				}
				branch.centerPos=Math.round((leaves.leaf1.y+leaves.leaf2.y)/2);
			}
			else { // column tree
				for (let leaf in leaves) {
					if (branch[leaf].heatCellList) { // label
						leaves[leaf].x=branch[leaf].aspect.attr('x');
						leaves[leaf].y=treeStart[type];
					}
					else { // branch
						branch[leaf].parent=branch;
						leaves[leaf].x=branch[leaf].centerPos;
						leaves[leaf].y=branch[leaf].connectionPos;
					}
				}
				branch.centerPos=Math.round((leaves.leaf1.x+leaves.leaf2.x)/2);
			}
			let branchFullSize=Math.round((branch.distance/treeLength[type])*treeSpace[type]),
				pathStrg='M'+(leaves.leaf1.x-0.5)+' '+(leaves.leaf1.y-0.5),
				bNode;
			if (type=='row') { // row tree
				branch.connectionPos=(rowLabelLocation=='left')? treeStart[type]+branchFullSize : treeStart[type]-branchFullSize;
				pathStrg+=' l'+(branch.connectionPos-leaves.leaf1.x)+' 0 l0 '+(leaves.leaf2.y-leaves.leaf1.y)+' l'+(leaves.leaf2.x-branch.connectionPos)+' 0';
				bNode=HM.paper.circle(branch.connectionPos-0.5,branch.centerPos-0.5,3);
			}
			else { // column tree
				branch.connectionPos=(colLabelLocation=='top')? treeStart[type]+branchFullSize : treeStart[type]-branchFullSize;
				pathStrg+=' l0 '+(branch.connectionPos-leaves.leaf1.y)+' l'+(leaves.leaf2.x-leaves.leaf1.x)+' 0 l0 '+(leaves.leaf2.y-branch.connectionPos);
				bNode=HM.paper.circle(branch.centerPos-0.5,branch.connectionPos-0.5,3);
			}
//console.log(pathStrg);
			//bNode.attr({stroke:colorList[branch.selectColorIdx],fill:colorList[branch.selectColorIdx]})
			bNode.attr({stroke:'none',fill:'#F00','fill-opacity':0})
			.data({branch:branch})
			//.hover(function(){this.attr({stroke:'#F00',fill:'#F00',cursor:'pointer'})},function(){var nodeColor=colorList[this.data('branch').selectColorIdx]; this.attr({stroke:nodeColor,fill:nodeColor,cursor:'auto'});})
			.hover(function(){this.transform('s1.7').attr({'fill-opacity':1,cursor:'pointer'});},function(){if(this.id) {this.attr({'fill-opacity':0,cursor:'auto'}).transform('');}}) //Tree is redrawn if exist sub node => node is deleted during hover!!! (problem for Chrome only so far)
			.click(function(e){if (e.ctrlKey || e.shiftKey) {rotateBranch(this.data('branch'));} else {highLightBranch(this.data('branch'));}});
			//.dblclick(function(){rotateBranch(this.data('branch'));})
			//.click(function(){highLightBranch(this.data('branch'));});
			branch.aspect.push(HM.paper.path(pathStrg).attr({stroke:colorList[branch.selectColorIdx],'stroke-width':1})); // [0] leaves
			branch.aspect.push(bNode); // [1] node
		}
		/* Sending all center circles to front */
		for (let i=0; i<branchList[type].length; i++) {branchList[type][i].aspect[1].toFront();}
		// Tree popups
		if (labelsOnSelect[type]) {
			let t=HM.paper.text(0.5,0.5,labelsOnSelect[type].text).attr({'font-weight':'bold','font-size':14,fill:'#FFF'}),
				b=t.getBBox(),
				tb=HM.paper.rect(b.x-2,b.y-2,b.width+4,b.height+5,2).attr({stroke:'#000',fill:'#000','stroke-opacity':0.7,'fill-opacity':0.7});
			t.toFront();
			treePopup[type]=HM.paper.set(tb,t).hide()
			.hover(function(){this.attr({cursor:'pointer'});},function(){this.attr({cursor:'auto'});})
			.click(function(){focussedLabelsAction('tree',type);}); // former treeClickAction
		}
	}

	function chooseColorIndex(usedColors) {
		let newColorIdx=null;
		for (let i=1; i<colorList.length; i++) {
			if (!usedColors[i]) {
				newColorIdx=i;
				usedColors[i]=true;
				break;
			}
		}
		if (newColorIdx===null) { // all colors used -> reset
			newColorIdx=1; // 1: already set to true
			//for (var i=2; i<colorList.length; i++) {delete usedColors[i];}
			usedColors={};
			usedColors[1]=true;
		}
		return newColorIdx;
	}

    function highLightBranch(branch,newColorIdx=null,usedColorIdx=null) {
		if (newColorIdx===null) {
			selectedLabels=[]; // reset list of selected labels
			if (branch.selectColorIdx > 0) { // unselect -> black
				newColorIdx=0;
				if (!branch.parent || branch.parent.selectColorIdx != branch.selectColorIdx) { // no parent or parent has another color / is not selected
					delete branchColorIndex[branch.type][branch.selectColorIdx]; // release color from used list
//console.log(branch.selectColorIdx,'released 1');
				}
				else {usedColorIdx=branch.selectColorIdx;}
				branchIsSelected=false;
			}
			else { // choose a new color
				newColorIdx=chooseColorIndex(branchColorIndex[branch.type]);
				branchIsSelected=true;
				//Show popup
				if (labelsOnSelect[branch.type]) {
					let translationStrg=(branch.type=='row')? 't'+branch.connectionPos+','+(branch.centerPos-15) : 't'+branch.centerPos+','+(branch.connectionPos-15);
					treePopup[branch.type][0].attr({stroke:colorList[newColorIdx],fill:colorList[newColorIdx]}); // box
					treePopup[branch.type].transform(translationStrg).data({trString:translationStrg}).toFront().show();
					setTimeout(function(){clearTreePopup(branch.type,translationStrg);},5000);
				}
			}
		}
		else if (branch.selectColorIdx && branch.selectColorIdx != usedColorIdx) {
			delete branchColorIndex[branch.type][branch.selectColorIdx]; // release color from used list (can have been already deleted by parent leaf)
//console.log(branch.selectColorIdx,'released 2');
		}
		let branchColor=colorList[newColorIdx];
		branch.aspect[0].attr({stroke:branchColor}); // leaves
		//branch.aspect[1].attr({fill:branchColor}); // node (Commented. Do not change color: always #F00)
		branch.selectColorIdx=newColorIdx;
		/* Propagating (de)selection */
		if (branch.leaf1.centerPos) { // points to another branch
			highLightBranch(branch.leaf1,newColorIdx,usedColorIdx);
		}
		else {
			highLightLabel(branch.leaf1,newColorIdx); // points to a label
			if (branchIsSelected) selectedLabels.push(branch.leaf1.id);
		}
		if (branch.leaf2.centerPos) {highLightBranch(branch.leaf2,newColorIdx,usedColorIdx);} // points to another branch
		else {
			highLightLabel(branch.leaf2,newColorIdx); // points to a label
			if (branchIsSelected) selectedLabels.push(branch.leaf2.id);
		}
    }

    function rotateBranch(branch) {
//console.log(branch.leaf1);
		let labelList={leaf1:[],leaf2:[]},numLabels={};
		//List of labels linked to leaf1
		for (let leaf in labelList) {
			if (branch[leaf].centerPos) {pushLabels(branch[leaf],labelList[leaf]);} // points to another branch
			else {labelList[leaf].push(branch[leaf]);} // points to a label
			numLabels[leaf]=labelList[leaf].length;
		}
//console.log(labelList);
		//Swapping leaves
		let cellSize,axis;
		if (branch.type=='row') {
			cellSize=hCellH;
			axis='y';
		}
		else {
			cellSize=hCellW;
			axis='x';
		}
		let shift={leaf1:numLabels.leaf2*cellSize,leaf2:-numLabels.leaf1*cellSize},
			shiftRanks={leaf1:numLabels.leaf2,leaf2:-numLabels.leaf1},
			movedLabels=[];
		for (let leaf in shift) {
			for (let i=0; i<labelList[leaf].length; i++) {
				let label=labelList[leaf][i],
					l=label.aspect;
				//l.attr(axis,l.attr(axis)+shift[leaf]);
				if (branch.type=='row') {l.attr({y:l.attr('y')+shift[leaf]});}
				else {l.attr({x:l.attr('x')+shift[leaf]});} //rotation moved to rank update
				//label.rank+=(leaf=='leaf1')? numLabels.leaf2 : -numLabels.leaf1;
 				movedLabels.push([label,shiftRanks[leaf]]);
				moveCells(label,axis,shift[leaf]);
			}
		}
		//Swapping labels & cells in rowList/columnList data structure
		let movedList,otherList,movedTypeIdx,otherTypeIdx,focusShift=0,focusShift2=0;
		if (branch.type=='row') {
			movedList=rowList;
			otherList=columnList;
			movedTypeIdx='rowIdx';
			otherTypeIdx='columnIdx';
			if (hasCells) {
				if (rowLabelLocation=='left') {
					focusShift=focusArea.attr('x')-cellX0;
					focusShift2=focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
				}
				else {
					focusShift=focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
					focusShift2=focusArea.attr('x')-cellX0;
				}
			}
		}
		else { // column
			movedList=columnList;
			otherList=rowList;
			movedTypeIdx='columnIdx';
			otherTypeIdx='rowIdx';
			if (hasCells) {
				if (colLabelLocation=='top') {
					focusShift=focusArea.attr('y')-cellY0;
					focusShift2=focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
				}
				else {
					focusShift=focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
					focusShift2=focusArea.attr('y')-cellY0;
				}
			}
		}
		let movedCells={},
			minLabelIdx=movedList.length-1,
			maxLabelIdx=0;
		for (let i=0; i<movedLabels.length; i++) {
			let label=movedLabels[i][0],
				rkShift=movedLabels[i][1],
				oldIdx=label.rank-1;
			label.rank+=rkShift;
			let newIdx=label.rank-1;
			for (let j=0; j<label.heatCellList.length; j++) { // update hCell (row/column)Idx
				label.heatCellList[j][movedTypeIdx]=newIdx;
			}
			movedList[newIdx]=label;
			movedCells[newIdx]=[];
			minLabelIdx=Math.min(minLabelIdx,newIdx);
			maxLabelIdx=Math.max(maxLabelIdx,newIdx);
			for (let j=0; j<otherList.length; j++) {
				movedCells[newIdx][j]=otherList[j].heatCellList[oldIdx];
			}

			/* Checking if label (& annotations) is inside focusArea */
			let trF='',trF2='';
			if (branch.type=='row') {
				if (hasCells && focusArea.data('status')=='on' && newIdx >= focusArea.data('labelIdx')[0] && newIdx <= focusArea.data('labelIdx')[1]) { // inside focus zone
					trF='t'+focusShift+',0';
					trF2='t'+focusShift2+',0';
					label.aspect.toFront(); // shift & move above bgRow
				}
				label.aspect.transform(trF); // also when empty in case moved outside focus
			}
			else { // column
				if (hasCells && focusArea.data('status')=='on' && newIdx >= focusArea.data('labelIdx')[2] && newIdx <= focusArea.data('labelIdx')[3]) { // inside focus zone
					trF='t0,'+focusShift;
					trF2='t0,'+focusShift2;
					label.aspect.toFront(); // move above bgColumn
				}
				label.aspect.transform(trF+'r'+colLabelAngle+','+label.aspect.attr('x')+','+label.aspect.attr('y'));
			}
			/* Moving annotation cells */
			for (let setName in label.annotations) {
				label.annotations[setName].aspect.toFront().transform(trF2);
			}
		}
		if (cellsValueTrace) {
			/* Redrawing flanking cells if column swap (not needed for row swap)*/
			if (branch.type != HM.normProcess.scope) {
				if (--minLabelIdx>=0) {moveCells(movedList[minLabelIdx],axis,0);} // redraw lower falnking cell to clean value lines
				else {minLabelIdx++;} // = 0
				if (++maxLabelIdx < movedList.length) {moveCells(movedList[maxLabelIdx],axis,0);} // redraw lower falnking cell to clean value lines
				else {maxLabelIdx--;} // = 0
			}
			/* Redrawing value lines */
			for (let i=minLabelIdx; i<=maxLabelIdx; i++) {
				drawCellsValueTrace(movedList[i].heatCellList);
			}
		}

		for (let idx in movedCells) {
			for (let i=0; i<movedCells[idx].length; i++) {
				otherList[i].heatCellList[idx]=movedCells[idx][i];
				otherList[i].heatCellList[idx][otherTypeIdx]=i;  // update hCell (otherType)Idx
			}
		}
		movedCells=null;

		//Swapping pointers as well
		let newLeaf2=branch.leaf1;
		branch.leaf1=branch.leaf2;
		branch.leaf2=newLeaf2;
//console.log(branch.leaf1);
		if (branch.leaf1.centerPos || branch.leaf2.centerPos) { // There a subBranches
			// Delete Tree
/*
			for (var i=0; i<branchList[branch.type].length; i++) {
				var delBrch=branchList[branch.type][i];
				delBrch.aspect[0].remove(); // path
				delBrch.aspect[1].unclick(); // circle
				delBrch.aspect[1].removeData();
				delBrch.aspect[1].remove();
			}
*/
			clearTree(branch.type);
			drawTree(branch.type);
		}
    }

    function pushLabels(branch,labelLeafList) {
		if (branch.leaf1.centerPos) {pushLabels(branch.leaf1,labelLeafList);} // points to another branch
		else {labelLeafList.push(branch.leaf1);} // points to a label

		if (branch.leaf2.centerPos) {pushLabels(branch.leaf2,labelLeafList);} // points to another branch
		else {labelLeafList.push(branch.leaf2);} // points to a label
    }


    /***********************************************************************/
    /***************************** HTML form *******************************/
    /***********************************************************************/
    function initializeHeatMapForm() {
		let htmlString='<TABLE cellspacing=0 cellpadding=0><TR>';
		if (HM.editMenu) {
			htmlString+='<TD nowrap><FIELDSET style="padding:2px"><LEGEND><B>Data normalization:</B></LEGEND>';
			/* Scope */
			if (HM.editMenu.scope) {
				htmlString+='<B>Scope:</B><SELECT onchange="idvis.registeredCharts['+HM.chartID+'].normProcess.scope=this.value; idvis.registeredCharts['+HM.chartID+'].computeHeatMap()"><OPTION value="row">'+HM.entities.row+'</OPTION><OPTION value="column"';
				if (HM.normProcess.scope=='column') {htmlString+=' selected';}
				htmlString+='>'+HM.entities.column+'</OPTION><OPTION value="all"';
				if (HM.normProcess.scope=='all') {htmlString+=' selected';}
				htmlString+='>Both</OPTION></SELECT>';
			}
			else if (HM.editMenu.scope==false && HM.normProcess.reference != 'z-score') {htmlString+='<B>Scope: '+HM.normProcess.scope+'</B>';}
			/* Reference */
			if (HM.editMenu.type) {
				htmlString+='&nbsp;&nbsp;<B>Reference:</B><SELECT onchange="idvis.registeredCharts['+HM.chartID+'].normProcess.reference=this.value; idvis.registeredCharts['+HM.chartID+'].computeHeatMap()"><OPTION value="z-score">z-score</OPTION><OPTION value="mean"';
				if (HM.normProcess.reference=='mean') {htmlString+=' selected';}
				htmlString+='>mean</OPTION><OPTION value="median"';
				if (HM.normProcess.reference=='median') {htmlString+=' selected';}
				htmlString+='>median</OPTION><OPTION value="mid-range"';
				if (HM.normProcess.reference=='mid-range') {htmlString+=' selected';}
				htmlString+='>mid-range</OPTION><OPTION value="min"';
				if (HM.normProcess.reference=='min') {htmlString+=' selected';}
				htmlString+='>min. value</OPTION><OPTION value="zero"';
				if (HM.normProcess.reference=='zero') {htmlString+=' selected';}
				var zeroStrg=(cellOnMouseOver)? cellOnMouseOver(0) : 'zero (0)';
				htmlString+='>'+zeroStrg+'</OPTION>';
				if (HM.normProcess.reference=='user') {
					htmlString+='<OPTION value="user" selected>user</OPTION>';
				}
				htmlString+='</SELECT>';
			}
			else if (HM.editMenu.type==false) { // nothing displayed if not declared at all
				htmlString+='&nbsp;&nbsp;<B>Reference: ';//+HM.normProcess.reference+'</B>';
				if (HM.normProcess.reference=='zero' && cellOnMouseOver) {htmlString+=cellOnMouseOver(0);}
				else {htmlString+=HM.normProcess.reference;}
				htmlString+='</B>';
			}
			/* Color */
			if (HM.editMenu.color) {
				htmlString+='&nbsp;&nbsp;<B>Color:</B><SELECT onchange="idvis.registeredCharts['+HM.chartID+'].normProcess.colors=this.value; idvis.registeredCharts['+HM.chartID+'].computeCellsColor()">';
				for (let i=0; i<colorPatterns.length; i++) {
					htmlString+='<OPTION value="'+colorPatterns[i][0]+'"';
					if (HM.normProcess.colors==colorPatterns[i][0]) {htmlString+=' selected';}
					htmlString+='>'+colorPatterns[i][1]+'</OPTION>';
				}
				htmlString+='</SELECT>';
			}
			/* Cell value trace */
			if (hCellW >= 20 || hCellH >= 20) {
				let checkStrg=(cellsValueTrace)? ' checked' : '';
				let disabStrg=((HM.normProcess.scope=='column' && hCellW < 20) || (HM.normProcess.scope != 'column' && hCellH < 20))? ' disabled' : '';
				htmlString+='&nbsp;&nbsp;<INPUT type="checkbox" id="'+formDivID+'_trace" '+checkStrg+' onchange="idvis.registeredCharts['+HM.chartID+'].showCellsValueLine(this.checked)"'+disabStrg+'/><B>Show value trace</B>&nbsp;';
			}
			htmlString+='</FIELDSET></TD>\n';
			/* Cell flag visibility */
			if (existFlaggedCells) {
				htmlString+='<TD nowrap><FIELDSET style="padding:2px"><LEGEND><B>Flag:</B></LEGEND><INPUT type="checkbox" checked onchange="idvis.registeredCharts['+HM.chartID+'].showCellFlag(this.checked)"/><B>Show '+flagText+'&nbsp;</FIELDSET></TD>\n';
			}
		}
/* Search */
if (HM.editMenu.searchable) {
    htmlString+='<FIELDSET style="padding:2px;white-space:nowrap"><LEGEND><B>Search:</B></LEGEND>';
    let searchBoxID=formDivID+'_search',
        searchResDivID=formDivID+'_srchRes';
    htmlString+='<INPUT type="text" id="'+searchBoxID+'" style="width:140px" value=""><INPUT type="button" value="Go" style="font-weight:bold;font-size:12px;width:40px" onclick="idvis.lib.searchDataPoints(idvis.registeredCharts['+mainC.chartID+'],document.getElementById(\''+searchBoxID+'\').value,\''+searchResDivID+'\',';
    if (typeof(HM.editMenu.searchable)=='object') {
        let extSearchText=HM.editMenu.searchable.text || 'Extended search',
        extSearchChkID=formDivID+'_extSearch';
        htmlString+='\''+extSearchChkID+'\')"><BR><INPUT type="checkbox" id="'+extSearchChkID+'" value="1">'+extSearchText;
    }
    else {
        htmlString+='null)">';
    }
    htmlString+='<DIV id="'+searchResDivID+'"></DIV>';
    htmlString+="</FIELDSET>\n";
}
		/* Export image button */
		if (HM.exportAsImage) {
			htmlString+='<TD>';
			if (HM.editMenu) htmlString+="&nbsp;&nbsp;";
            if (HM.exportAsImage[3]) { // image format
				htmlString+='<INPUT type="button" value="'+HM.exportAsImage[0]+'" style="font-weight:bold;font-size:10px" onclick="idvis.registeredCharts['+HM.chartID+'].mergeAndExportToImage(\''+HM.exportAsImage[3]+'\')"/></TD>\n';
			}
			else {
				htmlString+='<FONT style="font-weight:bold">'+HM.exportAsImage[0]+':</FONT>';
				//PNG
				htmlString+='<INPUT type="button" value="PNG" style="font-weight:bold;font-size:10px" onclick="idvis.registeredCharts['+HM.chartID+'].mergeAndExportToImage(\'png\')"/>';
				//SVG
				htmlString+='<INPUT type="button" value="SVG" style="font-weight:bold;font-size:10px" onclick="idvis.registeredCharts['+HM.chartID+'].mergeAndExportToImage(\'svg\')"/></TD>\n';
			}
        }
		//if (HM.editMenu) {htmlString+="</FIELDSET>\n";}
		htmlString+='</TR></TABLE>\n';

		document.getElementById(formDivID).innerHTML=htmlString;
    }

	/***********************************************************/
    /****************** Merge CANVAS with SGV ******************/
	/***********************************************************/
	this.mergeAndExportToImage=function(format)  {
		let tmpMap=this.paper.image(canvas.toDataURL(),cellX0,cellY0,mapW,mapH);
		tmpMap.toBack(); // So cell selection (if any) remains visible
		helpLogo.hide();
		helpPopup.hide();
		idvis.lib.exportSVGtoImg(paperDivID,this.exportAsImage[1],this.exportAsImage[2],format);
		tmpMap.remove();
		helpLogo.show();
		if (helpIsVisible) helpPopup.show();
	};

    /***********************************************************/
    /************************ Help text ************************/
    /***********************************************************/
    function showHideHelpText() {
		if (helpIsVisible) {
			helpPopup.animate({'fill-opacity':0,'stroke-opacity':0},300,'linear',function() {helpPopup.hide();});
			//helpPopup.hide();
			helpIsVisible=false;
		}
		else {
			helpPopup.show();
			let anim=Raphael.animation({'fill-opacity':0.7,'stroke-opacity':0.7},300,'linear');
			helpPopup[0].animate(anim); // box
			helpPopup[1].animateWith(helpPopup[0],anim,{'fill-opacity':1},300,'linear');

			helpIsVisible=true;
		}
    }



/*========================================= Nested objects ===========================================*/

	/*****************************************************************
					Heatmap cell object
	******************************************************************/
	function HeatCell(CH,value) {
//value*=-1;
		this.chart=CH;
		//this.value=(typeof(value)==='undefined' || value=='NA' || value==null)? null : value*1;
		this.value=(idvis.isNumber(value))? value : null; //isNaN('')=false!!!!
		this.pcValue=null;
		this.zscore=null;
		this.color=null;
		this.rowIdx=null;
		this.columnIdx=null;
		this.isExcluded=false;
		this.isFlagged=false;
		this.isSelectedData=null;
		this.aspect=null; // graphical object
	}

	/*****************************************************************
			Row/Column Label object
	******************************************************************/
	function Label(CH,type,rank,info) {
		this.chart=CH;
		this.type=type; //row or column
		this.rank=rank;
		this.label=info[0];
		this.id=(info[1])? info[1] : info[0];
		this.popupText=(info[2])? info[2] : null;
		this.heatCellList=[];
		this.annotations={};
		this.selectColorIdx=0;
		this.groupIndex=null;
		this.aspect=null; // graphical object
	}

	/*****************************************************************
			Group Label object
	******************************************************************/
	function labelGroup(CH,type,info) {
		this.chart=CH;
		this.type=type; //row or column
		this.label=info[0];
		this.id=(info[1])? info[1] : info[0];
		this.popupText=(info[2])? info[2] : null;
		this.labelIndexRange=[info[3],info[4]];
		//this.aspect; // graphical object
	}

	/*****************************************************************
					  Annotation set object
	******************************************************************/
	function AnnotationSet(type,name) {
		this.type=type; //row or column
		this.name=name;
		this.rank=null;
		this.popupText=name+'\nDoudle-click to delete';
		this.cells=[];
		this.aspect=null; // graphical object
	}

	/*****************************************************************
					  Annotation cell object
	******************************************************************/
	function AnnotationCell(set,name,label,color) {
		this.set=set;
		this.label=label;
		this.name=name;
		this.color=color; // not used
		this.aspect=null; // graphical object
	}

	/*****************************************************************
					  Branch object (Clustering tree piece)
	******************************************************************/
	function TreeBranch(CH,type,rank,leaf1,leaf2,distance) {
		this.chart=CH;
		this.type=type; //row or column
		this.rank=rank*1; // force number !!!rank of label at the end of all childen leaves1!!!
		this.leaf1=leaf1; // Label or Branch object
		this.leaf2=leaf2; // Label or Branch object
		this.distance=(distance)? distance : 1;
		this.selectColorIdx=0;
		// defined at draw -->
		this.parent=null;
		this.centerPos=null;
		this.connectionPos=null;
		this.aspect=[]; // array of graphical objects
	//console.log(this);
	}


};

/*
####>Revision history<####
# 3.2.0 [FEATURE] Support for rounded cells (PP 07/05/21)
# 3.1.6 [ENHANCEMENT] Improved annotation set color selection (PP 21/03/20)
# 3.1.5 [ENHANCEMENT] Improved handling of annotation values rank (PP 28/02/20)
# 3.1.4 Displays both PNG and SVG image export options if format is not specified by user (PP 06/06/18)
# 3.1.3 [Fix] Bug due to conflict between hCellsSensor.drag and .click in non-Firefox browsers (PP 09/02/18)
# 3.1.2 Added 2nd palette trim for single-color heatmap (PP 08/12:17)
# 3.1.1 Updates main DIV size when adding/removing annotations legends (PP 11/11/17)
# 3.1.0 Value trace management, labelOnMove() callback & new color patterns (PP 07/07/17)
# 3.0.0 Updated heatMap.js to code-encapsulated idvis.heatMap.js (PP 30/05/17)
# 2.3.3 Bug fix in translationStrg declaration during branch selection (PP 30/05/17)
# 2.3.2 Bug fixes in palette & column label emphasis for 1D clustering (11/01/17)
# 2.3.1 setGrgoups function replaced by defineGroups & other minor optimizations (PP 14/10/16)
# 2.3.0 Label group management (PP 12/10/16)
# 2.2.1 Improved label space computation using getBBox() & handles map/cell height/width parameters (PP 28/07/16)
# 2.2.0 Minimum heatCell height can be as low as 1px if HM.noMinCellHeight set to 'true' (PP 29/05/16)
# 2.1.9 Annotations cells also move around focus area (PP 28/05/16)
# 2.1.8 Uses chart registering for dynamic html calls (PP 20/02/16)
# 2.1.7 Max range & flag visibility options (PP 27/01/16)
# 2.1.6 Added check on label id match in addTree function (PP 15/01/16)
# 2.1.5 Flagged cells management (PP 10/10/15)
# 2.1.4 Bug fix for undefined hCell when re-dragging inside a focus area (PP 05/08/15)
# 2.1.3 Bug fix for mouse position on hCellsSensor using chartLibrary2 getMousePositionInElement() (PP 24/07/15)
# 2.1.2 Bug fixes in label swap in 1D view & multiple values per annotation cell (PP 28/04/15)
# 2.1.1 Bug fix in focus-based label list extraction (PP 20/03/15)
# 2.1.0 Drag focus area & replacement of treeOnClick with labelsOnSelect (PP 10/03/15)
# 2.0.3 Loose id match option in addAnnotationCell() & fixes bugs in display and in tree color release (PP 09/01/15)
# 2.0.2 Palette with values distribution (PP 21/12/14)
# 2.0.1 Compatibility with 1D-vertical clustering (PP 27/11/14)
# 2.0.0 Uses CANVAS for heatmap cells. Requires chartLibrary2.js only if exportAsImage is defined (PP 03/10/14)
# 1.5.4 Returns a success status on addAnnotation() & optional custom color argument: 'comma-separated matchedList[::#color]' (PP 05/08/14)
# 1.5.3 Added option for exporting set of selected row/column label ids to callback function (PP 08/07/14)
# 1.5.2 Popup on annotation labels, hidden tree nodes & other minor improvements (PP 08/07/14)
# 1.5.1 Calls annotation removal callback function (if any) before removing annotation (PP 17/06/14)
# 1.5.0 Annotation removal & new heatCell event management (PP 14/06/14)
# 1.4.2 cell update on label-level (ex/in)clusion (PP 17/04/14)
# 1.4.1 noCellExclusion flag & no-value cells can be excluded (PP 11/04/14)
# 1.4.0 Transmission of cell exclusion event & pre-exclusion at start (PP 10/04/2014)
# 1.3.9 Added getDivID function && optional export image button (PP 11/02/14)
# 1.3.8 GPL license (PP 23/09/13)
# 1.3.7 this.renameLabel function & allow no-value cells (de)selection on row/column double click events (PP 13/05/13)
# 1.3.6 No-value cells not selected on row/column double click events (PP 26/03/13)
# 1.3.5 Added dblclick & shift+dblclick events on row/column labels for multiple cell selection (PP 25/03/13)
# 1.3.4 Added clearAllCells function (PP 19/03/13)
# 1.3.3 Minor bug fix in addAnnotation function (PP 29/01/13)
# 1.3.2 Bug fix in function updateCellsValue (PP 10/12/12)
# 1.3.1 Rank of Tree branch based on label rank: no longer uses addTree() data[0] (PP 08/12/12)
# 1.3.0 ModClick on labels & bug fix (PP 06/12/12)
# 1.2.9 modKeyPressed declared by chartLibrary.js (PP 03/12/12)
# 1.2.8 Handles user-defined limit values: min,ref,max (PP 28/11/12)
# 1.2.7 Minor change in annotation popup (PP 02/10/12)
# 1.2.6 Added annotation stripes (PP 11/09/12)
# 1.2.4 Bug fix for cell hover handled by external function (PP 14/06/12)
# 1.2.3 Branch parent relationship & fix minor bug in branch color selectection (PP 06/06/12)
# 1.2.2 Multiple color selection for tree branches (PP 05/06/12)
# 1.2.1 Help string for clustering tree (PP 04/06/12)
# 1.2.0 Added clustering tree 1D or 2D (PP 04/06/12)
# 1.1.1 Help popup text (PP 26/05/12)
# 1.1.0 Cell middle X & bottom Y are sent upon selection (PP 22/05/12)
# 1.0.9 More z-score debuging (PP 14/05/12)
# 1.0.8 z-score correction (PP 12/05/12)
# 1.0.7 Shift + mouseover for more details on cell (PP 10/05/12)
# 1.0.6 Fixed bugs & more checks (PP 09/05/12)
# 1.0.5 Added form menu (PP 06/05/12)
# 1.0.4 Multiple color ranges & multi-cell selection (PP 05/05/12)
# 1.0.3 draggable labels (PP 29/04/12)
*/
