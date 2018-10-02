/*
################################################################################
# heatMap.js        2.3.5                                                      #
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
                Heatmap cell object
******************************************************************/
function HeatCell(CH,value) {
    this.chart=CH;
    //this.value=(typeof(value)==='undefined' || value=='NA' || value==null)? null : value*1;
	this.value=(value==null || value.length==0 || isNaN(value))? null : value; //isNaN('')=false!!!!
    this.pcValue;
    this.zscore;
    this.color;
    //this.row;
    //this.column;
	this.rowIdx;
	this.columnIdx;
    this.isExcluded=false;
	this.isFlagged=false;
    this.isSelectedData;
    this.aspect; // graphical object
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
    this.aspect; // graphical object
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
	this.rank;
	this.popupText=name;
	this.cells=[];
	this.aspect; // graphical object
}
/*****************************************************************
                  Annotation cell object
******************************************************************/
function AnnotationCell(set,name,label,color) {
	this.set=set;
	this.label=label;
	this.name=name;
	this.color=color; // not used
    this.aspect; // graphical object
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
	this.parent;
    this.centerPos;
    this.connectionPos;
    this.aspect=[]; // array of graphical objects
//console.log(this);
}


/*****************************************************************
                  HeatMap object
******************************************************************/
function heatMap(mapData) {
    var HM=this;
	this.chartID=registerInCuBioJS(this);

    /* Configuration parameters */
    var mainDivID=mapData.div; //document.getElementById(mapData.div);
    var mainDiv=document.getElementById(mainDivID);
	this.exportAsImage=mapData.exportAsImage;
    var paperDivID=mainDivID+'_paper';
    var canvasID=mainDivID+'_canvas';
	var formDiv;
	this.editMenu;
    if (mapData.editable || mapData.editableItems || this.exportAsImage) { // null if not declared
		this.editMenu=(mapData.editableItems)? mapData.editableItems : (mapData.editable)? {scope:true, type:false, color:true} : null; // used by external JS to know is form exists
		formDivID=mainDivID+'_form';
    }
    else {
		this.editMenu=null;
    }
	this.getDivID=function(){return paperDivID;}
    this.normProcess={ // default
		scope:'row', // 'row' 'column' 'all'
		reference:'mean', // z-score,median,mean,mid-range,min,max,user
		colors:'GBR', // GBR,GYR,GWR,BBR,BYR,BWR
		colorList:{limit:{'G':'#0F0','R':'#F00','B':'#00F'},reference:{'B':'#000','Y':'#FF0','W':'#FFF'}}
    };
	var maxRange=mapData.maxRange || {}; // max range for cell values {type:<absolute|quartile>,value:<[min,max]|number of quartiles>}
	var rangeMinValue,rangeMaxValue;
	if (mapData.entities) { // Which entities are reprensented by row/column
		this.entities=mapData.entities;
		if (!this.entities.row) this.entities.row='Row';
		if (!this.entities.column) this.entities.row='Column';
	}
	else {this.entities={row:'Row',column:'Column'};}
    if (mapData.normalization) {
		for (var attr in mapData.normalization) {this.normProcess[attr] = mapData.normalization[attr];}
		if (!this.normProcess.scope.match('row|column|all')) {this.normProcess.scope='row';}
		if (!this.normProcess.reference.match('mean|median|mid-range|min|zero|user') || (this.normProcess.reference=='user' && !this.normProcess.limitValues)) {this.normProcess.reference='z-score';}
		if (!this.normProcess.colors.match('GBR|GYR|GWR|BBR|BYR|BWR')) {this.normProcess.colors='GBR';}
		if (this.normProcess.reference=='z-score') {this.editMenu.scope=false;}
		//else if (this.normProcess.reference=='user') {this.editMenu.type=false;}
    }
//this.treeSpace={row:0,column:0};
    var rowLabelLocation=mapData.rowLabelPos || 'left';
    var colLabelLocation=mapData.colLabelPos || 'top';
	var colLabelOrientation=mapData.colLabelOrientation || 'vertical';
    var colLabelAngle=(colLabelOrientation=='diagonal')? -45 : -90;
	var colLabelRatio=Math.abs(Math.sin(colLabelAngle * (Math.PI / 180))); //radian
	var movableLabel={};
	var selectedLabels=[]; // list of label IDs highlighted by clicking on a tree node
	var minCellH=(mapData.noMinCellHeight)? 1 : 5;
    movableLabel.row=(mapData.moveLabel || mapData.moveRow)? true : false;
    movableLabel.column=(mapData.moveLabel || mapData.moveColumn)? true : false;

	//Palette
 	var palette={valueDistrib:[],trimValues:{min:100,max:100,svgMin:null,svgMax:null},hasNegValues:false,hasPosValues:false,posX:null,posY:null,refLimit:null,svgElements:[]};
	//palette.posX=(rowLabelLocation=='left')?

	// Row,column mouse events
    var rowOnClick=mapData.rowOnClick || null;
	var rowOnModClick=mapData.rowOnModClick || null;
	var rowOnDblClick=mapData.rowOnDblClick || null;
	var rowOnModDblClick=mapData.rowOnModDblClick || null;
	var rowGroupOnClick=mapData.rowGroupOnClick || null;
    var columnOnClick=mapData.columnOnClick || null;
	var columnOnModClick=mapData.columnOnModClick || null;
	var columnOnDblClick=mapData.columnOnDblClick || null;
	var columnOnModDblClick=mapData.columnOnModDblClick || null;
	var columnGroupOnClick=mapData.columnGroupOnClick || null;
	// Cell mouse events
	var oneCellSelection=(mapData.singleCellSelection)? true : false;
	var lastSelectedCell,lastHighlightedCell; // js objects
    var cellOnClick=mapData.cellOnClick || null;
	var noCellExclusion=(mapData.noCellExclusion)? true : false;
    var cellOnModClick=(!noCellExclusion && mapData.cellOnModClick)? mapData.cellOnModClick : null; // noCellExclusion takes over cellOnModClick
    var cellOnMouseOver=mapData.cellOnMouseOver || null;
	var annotationOnDelete=mapData.annotationOnDelete || null;
    var cellHighlightBgColor='#000'; //(this.normProcess.colors.match('Y'))? '#FFF' : '#000';
    var cellHighlightColor,noDataCellColor; // computed by computeCellsColor
	// Tree mouse event
	var labelsOnSelect=mapData.labelsOnSelect || {}; // {row/column:{text:...,action:...[,ranking:true/false]}}
	for (var type in {row:1,column:1}) {
		if (labelsOnSelect[type]) {
			if (!labelsOnSelect[type].action) {delete labelsOnSelect[type];}
			else if (!labelsOnSelect[type].text) {labelsOnSelect[type].text='Click for action';}
		}
	}
	// Flag
	var flagText=mapData.flagText || 'Flagged';
	var cellFlagVisible=true;
	var existFlaggedCells=false; // default. update in addRow()

	var treePopup={}; // row/column popup SVG
	var branchIsSelected=false; // Flag for tree branch selection
	var colorList=['#000000','#0000FF','#4AA02C','#F660AB','#FBB917','#8AFB17','#9A9A9A','#7E2217','#95B9C7','#E18B6B']; // for clustering tree/label/annotation selection
	var branchColorIndex={row:{},column:{}};
	//var annotColorIndex={row:{},column:{}};
	var annotNameToColor={row:{},column:{}}

    /*** Variables ***/
    /* Data */
    var columnList=[],rowList=[],groupList={row:[],column:[]};
	this.getColumnList=function() {return columnList;} // to allow access from outside HM
	this.getRowList=function() {return rowList;} // to allow access from outside HM
	var hasCells=false; // 1D or 2D
    var branchList={}; // Clustering tree {row:[],column:[]}
    var treeLength={};
	var annotationList={row:{},column:{}};
	var annotationCellSize={row:null,column:null};
	var annotationLegends; // SVG set containing all SVG elements of legends
    //var numDefinedValues=0;
	var maxColLabel='',maxRowLabel='',maxGroupAndLabels={column:'',row:''};
	//var maxColLabelSize=15; //7; // min.length used
	//var maxRowLabelSize=15; //7; // min.length used
   var maxAbsZscore=0;
    /* Graphics */
	var topSpace=20;
	var bottomSpace=20;
	var leftSpace=20;
	var rightSpace=20;
	var paperWidth,paperHeight,hCellW,hCellH,flagSize,rowLabelFontSize,colLabelFontSize,rowLabelSpace,colLabelSpace,colDiagSpace,mapW,mapH,cellX0,cellY0,canvas,mapContext;
	var paperWidth0,paperHeight0,background,backgroundPath;
	var treeSpace={row:50,column:50};
	var treeStart={row:0,column:0};
    //var mapGeometry={};
    var labelPopupText,cellPopupText;
	var hCellsSensor; //SVG rect over heatCells
	var focusArea; //SVG drag rect over heatCells
    var labelPosRanges={};
    var maxOverlap={};
    var dragParams={};
    var dragContext=false;
    var dblClickContext=false;
    var helpIsVisible=false;
    var helpLogo,helpPopup;

	/********** Set all columns ***************/
	this.setColumns=function(colData) {
		var rank=0;
		for (var i=0; i<colData.length; i++) {
			var newCol=new Label(this,'column',++rank,colData[i])
			columnList.push(newCol);
			if (maxColLabel.length < newCol.label.length) maxColLabel=newCol.label;
		}
//console.log('COL: '+columnList.length);
	}

	/********** Add multiple or single row(s) & associated data ***************/
	this.addRows=function(rowsData) {
		for (var i=0; i<rowsData.length; i++) {
			HM.addRow(rowsData[i][0],rowsData[i][1],rowsData[i][2],rowsData[i][3]);
		}
	}
	this.addRow=function(rowInfo,rowData,excludedCells,flaggedCells) {
		var newRow=new Label(this,'row',rowList.length+1,rowInfo);
		rowList.push(newRow);
		if (maxRowLabel.length < newRow.label.length) maxRowLabel=newRow.label;
//console.log('ROW: '+rowData.length);
		for (var i=0; i<columnList.length; i++) { // columnList.length in case less values in current row
			var hCell=new HeatCell(this,rowData[i]);
			//hCell.row=newRow;
			hCell.rowIdx=rowList.length-1;
			//hCell.column=columnList[i];
			hCell.columnIdx=i;
			columnList[i].heatCellList.push(hCell);
			newRow.heatCellList.push(hCell);
			if (hCell.value != null) {
				rangeMinValue=(rangeMinValue==null)? hCell.value : Math.min(rangeMinValue,hCell.value);
				rangeMaxValue=(rangeMaxValue==null)? hCell.value : Math.max(rangeMaxValue,hCell.value);
			}
		}
		// Excluded cells
		if (excludedCells) {
			for (var e=0; e<excludedCells.length; e++) {
				newRow.heatCellList[excludedCells[e]].isExcluded=true;
			}
		}
		// Flagged cells
		if (flaggedCells) {
			for (var e=0; e<flaggedCells.length; e++) {
				newRow.heatCellList[flaggedCells[e]].isFlagged=true;
				existFlaggedCells=true;
			}
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
		var labelList=(type=='row')? rowList : columnList;
		for (var i=0; i<groupData.length; i++) {
			var groupIdx=groupList[type].length;
			var gr=new labelGroup(this,type,groupData[i]);
			groupList[type].push(gr);
			var maxLabelIdx=gr.labelIndexRange[0];
			for (var l=gr.labelIndexRange[0]; l<=gr.labelIndexRange[1]; l++) {
				labelList[l].groupIndex=groupIdx;
				if (labelList[maxLabelIdx].label.length < labelList[l].label.length) maxLabelIdx=l;
			}
			gr.labelIndexRange.push(maxLabelIdx);
			var maxLabel=gr.label+labelList[maxLabelIdx].label;
//console.log(maxLabelIdx,maxLabel);
			if (maxGroupAndLabels[type].length < maxLabel.length) maxGroupAndLabels[type]=maxLabel;
		}
	}


	/********** Add a row or column clustering tree ***************/
	this.addTree=function(type,dataStrg) {
		movableLabel[type]=false; // overwrites (just to be safe)
	    branchList[type]=[];
	    treeLength[type]=0;
	    var labelList=(type=='row')? rowList : columnList;
	    var labelIdList={}; // id -> Label
	    for (var i=0; i<labelList.length; i++) {
			labelIdList[labelList[i].id]=labelList[i];
	    }
	    var treeData=dataStrg.split(';');
	    for (var i=0; i<treeData.length; i++) {
			var branchData=treeData[i].split(',');
			//var leaf1=(branchData[1]>0)? branchList[type][branchData[1]-1] : labelIdList[-branchData[1]];
			//var leaf2=(branchData[2]>0)? branchList[type][branchData[2]-1] : labelIdList[-branchData[2]];
			var leaf1=(branchData[1].match('-'))? labelIdList[branchData[1].substring(1)] : branchList[type][branchData[1]-1];
			var leaf2=(branchData[2].match('-'))? labelIdList[branchData[2].substring(1)] : branchList[type][branchData[2]-1];
			//Checks if label ids are matched
			if (!leaf1 && branchData[1].match('-')) {alert('ERROR: id "'+branchData[1].substring(1)+'" not found in '+type+' labels while reading tree.'); break;}
			if (!leaf2 && branchData[2].match('-')) {alert('ERROR: id "'+branchData[2].substring(1)+'" not found in '+type+' labels while reading tree.'); break;}
			if (leaf1.rank > leaf2.rank) { // invert leaves
				var tmpLeaf=leaf1;
				leaf1=leaf2;
				leaf2=tmpLeaf;
			}
			//branchList[type].push(new TreeBranch(this,type,branchData[0],leaf1,leaf2,branchData[3]));
			branchList[type].push(new TreeBranch(this,type,leaf1.rank,leaf1,leaf2,branchData[3]));
			treeLength[type]=Math.max(treeLength[type],branchData[3]);
	    }
	}

	/********************* Heat map display (Raphael) *********************/
    this.draw=function() {
		if (branchList.row) {treeSpace.row=Math.min(500,Math.max(150,Math.round(5*branchList.row[branchList.row.length-1].distance/branchList.row[0].distance)));} /* [100 >= 5px for minDist <= 500] */
		if (branchList.column) {treeSpace.column=Math.min(500,Math.max(150,Math.round(5*branchList.column[branchList.column.length-1].distance/branchList.column[0].distance)));} /* [100 >= 5px for minDist <= 500] */
//console.log(treeSpace);
		var maxMapWidth=(treeSpace.row)? 1500-treeSpace.row : 1500;
		//Columns
		var numColumns=columnList.length;
		//hCellW=Math.round(maxMapWidth/numColumns);
		//if (hCellW < 11) {hCellW=11;} else if (hCellW > 50) {hCellW=50;}
		hCellW=(mapData.cellWidth)? mapData.cellWidth : (mapData.width)? mapData.width/numColumns : Math.max(11,Math.min(50,Math.round(maxMapWidth/numColumns))); // mapData.cellWidth overwrites mapData.width!!!
		maxOverlap.column=hCellW/2;
		//Rows
		var numRows=rowList.length;
		if (numRows) { // 2D
			hasCells=true;
			//hCellH=Math.round(750/numRows);
			//if (hCellH < minCellH) {hCellH=minCellH;} else if (hCellH > 25) {hCellH=25;} // 13
			//if (hCellH > hCellW) {hCellH=hCellW;}
			hCellH=(mapData.cellHeight)? mapData.cellHeight : (mapData.height)? mapData.height/numRows : Math.min(hCellW,Math.max(minCellH,Math.min(25,Math.round(750/numRows))));
//hCellH=1;
//console.log(hCellH);
			maxOverlap.row=hCellH/2;
		}
		else {hCellH=0;} // 1D
		flagSize=(hCellH >= 5)? Math.min(10,Math.min(hCellW,hCellH)) : 3; //Size of flag triangle
//console.log('FLAG',flagSize);
//var rowTreeSpace=(branchList.row)? treeSpace : 0;
//var columnTreeSpace=(branchList.column)? treeSpace : 0;

		//Column labels space
		colLabelFontSize=(hCellW >= 25)? 14 : (hCellW >= 20)? 12 : 8;
		//Row labels space
		rowLabelFontSize=(hCellH >= 25)? 14 : (hCellH >= 20 || !hasCells)? 12 : 8; // 12 default size if 1D

		/* Temporary paper to compute label space */
		var tmpPaper=Raphael(mainDivID,0,0);
		var colText=(maxGroupAndLabels.column && maxGroupAndLabels.column.length > maxColLabel.length)? maxGroupAndLabels.column : maxColLabel;
		var colLabel=tmpPaper.text(0,0,colText).attr({'font-size':colLabelFontSize,'font-weight':'bold'});
		colLabelSpace=colLabel.getBBox().width;
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
			var rowText=(maxGroupAndLabels.row && maxGroupAndLabels.row.length > maxRowLabel.length)? maxGroupAndLabels.row : maxRowLabel;
			var rowLabel=tmpPaper.text(0,0,rowText).attr({'font-size':rowLabelFontSize,'font-weight':'bold'});
			rowLabelSpace=Math.max(rowLabel.getBBox().width,200); // room for palette
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
		var mainDiv=document.getElementById(mainDivID);
		var formDivHeight=0;
		if ((this.editMenu && hasCells) || this.exportAsImage) { // hasCells to check for 2D
			mainDiv.innerHTML='<DIV id="'+formDivID+'" style="margin-bottom:4px"></DIV>';
			/* Form menu */
			initializeHeatMapForm();
			formDivHeight=document.getElementById(formDivID).offsetHeight;
		}
		var divStrg='<DIV style="position:relative;">';
		divStrg+='<CANVAS id="'+canvasID+'" style="position:absolute;top:'+cellY0+'px;left:'+cellX0+'px" width="'+mapW+'" height="'+mapH+'"></CANVAS><DIV id="'+paperDivID+'" style="background:rgba(255,255,255,0);"></DIV>'; // class="transparent"
//divStrg+='<CANVAS id="'+canvasID+'_zoom" style="position:absolute;background-color:#0000FF;top:'+cellY0+'px;left:'+cellX0+'px" width="'+(mapW/2)+'" height="'+mapH+'"></CANVAS><DIV id="'+paperDivID+'_zoom" style="background:rgba(255,255,255,0);background-color:#00FFFF">Hello Patrick</DIV>'; // class="transparent"
		divStrg+='</DIV>';
		mainDiv.innerHTML+=divStrg;

		var paperDiv=document.getElementById(paperDivID);
		paperDiv.style.position='absolute'; // inside parent DIV (because relative)
		paperDiv.style.top='0px', paperDiv.style.left='0px';

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
		var backgroundPathStrg='M0,0 L'+paperWidth+',0 L'+paperWidth+','+paperHeight+' L0,'+paperHeight+' L0,'+cellY0+' L'+cellX0+','+cellY0+' L'+cellX0+','+(cellY0+mapH)+' L'+(cellX0+mapW)+','+(cellY0+mapH)+' L'+(cellX0+mapW)+','+cellY0+' L0,'+cellY0+' Z';
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

			/* Sensor */
			hCellsSensor=HM.paper.rect(cellX0,cellY0,mapW,mapH,0).attr({fill:'#555','fill-opacity':0,stroke:'none'})
			.mousemove(function(e){
							var hCell=findHeatCell(e);
							//if (!hCell.aspect) return; // mouseout too fast **********************************
							if (lastHighlightedCell) {
								if (hCell.rowIdx==lastHighlightedCell.rowIdx && hCell.columnIdx==lastHighlightedCell.columnIdx) {return;} //same cell
								else {setCellEmphasis(lastHighlightedCell,'off');} // new cell
							}
							setCellEmphasis(hCell,'on',e);
							lastHighlightedCell=hCell;
						})
			.mouseout(function(){if (lastHighlightedCell) {setCellEmphasis(lastHighlightedCell,'off'); lastHighlightedCell=null;}})
			.drag(function(dx,dy,x,y,e){dragFocusMove(dx,dy,x,y,e);}, //_extendDragging,
				  function(x,y,e){dragFocusStart(x,y,e);}, //_setDragging,
				  function(){dragFocusEnd();} //_endDragging
				  );
			if (cellOnClick) {hCellsSensor.click(function(e){if ((e.ctrlKey || e.shiftKey) && !noCellExclusion) {excludeCell(lastHighlightedCell);} else if (lastHighlightedCell) {selectCell(lastHighlightedCell);}});}
			else if (!noCellExclusion) {hCellsSensor.click(function(e){if (e.ctrlKey || e.shiftKey) {excludeCell(lastHighlightedCell);}});}

			/* Focus area & control(s) */
			var setY=7;
			var closeSet=HM.paper.set( // Close icon
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
				var selRowsSet=HM.paper.set( // Row list icon
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
				var selColumnsSet=HM.paper.set( // Column list icon
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
		var helpY=(colLabelLocation=='top')? paperHeight-12 : 12;
		var helpLink=this.paper.text(6,helpY,'HELP').attr({'text-anchor':'start','font-size':12,'font-weight':'bold',fill:'#FFF'});
		//helpLink.click(showHideHelpText);
		var hlBox=helpLink.getBBox();
		var hlb=this.paper.rect(hlBox.x-2,hlBox.y-1,hlBox.width+4,hlBox.height+2,0).attr({stroke:'none',fill:'#444','fill-opacity':1});
		//hlb.click(showHideHelpText);
		helpLogo=this.paper.set(helpLink,hlb);
		helpLink.toFront();
		this.paper.set(helpLink,hlb)
		.hover(function(){this.attr({cursor:'pointer'});},function(){this.attr({cursor:'auto'});})
		.click(showHideHelpText);
		var yShift=50;
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
		var helpTextY=(colLabelLocation=='top')? paperHeight-yShift : yShift;
		var hText=this.paper.text(10,helpTextY,helpText).attr({'text-anchor':'start','font-size':12,fill:'#FFF'});
		var htBox=hText.getBBox();
		var htb=this.paper.rect(htBox.x-6,htBox.y-4,htBox.width+12,htBox.height+8,10).attr({stroke:'none',fill:'#000'}); //,'fill-opacity':0.8
		hText.toFront();
		helpPopup=this.paper.set(htb,hText).attr({'fill-opacity':0,'stroke-opacity':0}).hide().click(showHideHelpText);

		/*** Drawing labels ***/
		function drawLabels() {
//console.time('DRAW_LABELS');
		    /** Rows **/
		    var rowX,rowLabelAnchor;
		    if (rowLabelLocation=='left') {
				rowX=leftSpace+rowLabelSpace-5;
				rowLabelAnchor='end';
		    }
		    else { // right
				rowX=leftSpace+treeSpace.row+mapW+5;
				rowLabelAnchor='start';
		    }
		    if (colLabelLocation=='bottom' && colLabelOrientation=='diagonal') {rowX+=colDiagSpace;}
		    var rowY=topSpace+Math.round(hCellH/2);
		    if (colLabelLocation=='top') {rowY+=colLabelSpace;}
		    else {rowY+=treeSpace.column;} // 0 or 200
			for (var i=0; i<rowList.length; i++) {
				var label;
				if (hCellH >= 5) {
					label=HM.paper.text(rowX,rowY,rowList[i].label)
					.attr({'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':rowLabelAnchor,'fill':colorList[rowList[i].selectColorIdx]})
					.data('row',rowList[i])
					.hover(function(){setLabelEmphasis(this,'on')},function(){setLabelEmphasis(this,'off')}) //;
					.click(function(e){if (dragContext) {dragContext=null;} else {var l=this; setTimeout(function(){labelClick(l,e)},300)}}) //setTimeout(labelClick,300,this,e) not compatible with IE9
					.dblclick(function(e){labelDoubleClick(this,e)});
					if (movableLabel.row) {label.drag(dragLabelMove,dragLabelStart,dragLabelStop);}
				}
				else { // no visible labels
					label=HM.paper.text(rowX,rowY,'').data('row',rowList[i]);
				}
				rowList[i].aspect=label;
				if (i==0) {labelPosRanges.minY=label.attr('y');}
				else if (i==rowList.length-1) {labelPosRanges.maxY=label.attr('y');}
				rowY+=hCellH;
		    }
			/* Row label groups */
			var shX=(rowLabelLocation=='left')? 5 : -5;
			var shY=hCellH/2;
			for (var i=0; i<groupList.row.length; i++) {
				var gr=groupList.row[i];
				var startY=rowList[gr.labelIndexRange[0]].aspect.attr('y');
				var endY=rowList[gr.labelIndexRange[1]].aspect.attr('y')
				var b=rowList[gr.labelIndexRange[2]].aspect.getBBox(); //longest label
				var barX,grX;
				if (rowLabelLocation=='left') {
					barX=b.x - 7.5;
					grX=barX - 7.5;
				}
				else {
					barX=b.x + b.width + 7.5;
					grX=barX + 7.5;
				}
				var bar=HM.paper.path('M'+(barX+shX)+','+(startY-shY)+' L'+barX+','+startY+' L'+barX+','+endY+' L'+(barX+shX)+','+(endY+shY)).attr({stroke:'#000','stroke-width':1});
				var grLabel=HM.paper.text(grX,startY+((endY-startY)/2),gr.label)
					.attr({'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':rowLabelAnchor,'fill':'#000'})
					.data({group:gr,bar:bar})
					.hover(function(){setLabelEmphasis(this,'on')},function(){setLabelEmphasis(this,'off')});
				if (rowGroupOnClick) {grLabel.click(function(e){rowGroupOnClick(this.data('group').id)});}
			}

//console.log('minY='+labelPosRanges.minY+', maxY='+labelPosRanges.maxY);
		    /** Columns **/
		    var colX=leftSpace+Math.round(hCellW/2);
		    if (rowLabelLocation=='left') {colX+=rowLabelSpace;}
		    else {colX+=treeSpace.row;} // 0 or 200
		    var colY,colLabelAnchor;
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
		    for (var i=0; i<columnList.length; i++) { // columns
				var label=HM.paper.text(colX,colY,columnList[i].label).
				attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':colLabelAnchor,'fill':colorList[columnList[i].selectColorIdx]})
				//.rotate(colLabelAngle,colX,colY)
				.transform('r'+colLabelAngle+','+colX+','+colY)
				.data('column',columnList[i])
				.hover(function(){setLabelEmphasis(this,'on')},function(){setLabelEmphasis(this,'off')})
				.click(function(e){if (dragContext) {dragContext=null;} else {var l=this; setTimeout(function(){labelClick(l,e)},300)}}) //setTimeout(labelClick,300,this,e) not compatible with IE9
				.dblclick(function(e){labelDoubleClick(this,e)});
				if (movableLabel.column) label.drag(dragLabelMove,dragLabelStart,dragLabelStop);
				columnList[i].aspect=label;
//if (i==0) console.log('DRAW: '+columnList[i].aspect.transform(),' <- ','r'+colLabelAngle+','+colX+','+colY);  //.toString()
				if (i==0) {labelPosRanges.minX=label.attr('x')} // rotation x->y
				else if (i==columnList.length-1) {labelPosRanges.maxX=label.attr('x')}
				colX+=hCellW;
		    }
			/* Column label groups */
			var shX=hCellW/2,shX2,shY;
			if (colLabelLocation=='top') {
				shX2=(colLabelOrientation=='diagonal')? 7.5 : 0;
				shY=5;
			}
			else {
				shX2=(colLabelOrientation=='diagonal')? -7.5 : 0;
				shY=-5;
			}
//console.log(colLabelFontSize,columnList[groupList.column[0].labelIndexRange[0]].aspect.attr('font-size'));
			for (var i=0; i<groupList.column.length; i++) {
				var gr=groupList.column[i];
				var b=columnList[gr.labelIndexRange[2]].aspect.getBBox(); //longest label
//HM.paper.rect(b.x,b.y,b.width,b.height);
				var startX,endX;
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
				var barY,grY;
				if (colLabelLocation=='top') {
					barY=b.y - 7.5;
					grY=barY - 7.5;
				}
				else {
					barY=b.y2 + 7.5;
					grY=barY + 7.5;
				}
				var bar=HM.paper.path('M'+(startX-shX-shX2)+','+(barY+shY)+'L'+startX+','+barY+' L'+endX+','+barY+' L'+(endX+shX-shX2)+','+(barY+shY)).attr({stroke:'#000','stroke-width':1});
				var grX=startX+(shX2)+(endX-startX)/2;
				var grLabel=HM.paper.text(grX,grY,gr.label)
					.attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':colLabelAnchor,'fill':'#000'})
					.transform('r'+colLabelAngle+','+grX+','+grY)
					.data({group:gr,bar:bar})
					.hover(function(){setLabelEmphasis(this,'on')},function(){setLabelEmphasis(this,'off')});
				if (columnGroupOnClick) {grLabel.click(function(e){columnGroupOnClick(this.data('group').id)});}
			}

		}
    }


	/*************************************************************/
    /********************** Computing heatmap ********************/
    /*************************************************************/
	function applyMaxRange() {
		function sortNumber(a,b) {return a-b;}
		function computeQuartile(ordValues,pos) {
			var calcValue;
			if (pos==Math.floor(pos)) {calcValue=ordValues[pos-1];}
			else {
				pos=Math.floor(pos);
				calcValue=(ordValues[pos-1] + ordValues[pos]) / 2;
			}
			return calcValue;
		}
		if (maxRange.type.match(/^abs/)) { // absolute
			rangeMinValue=maxRange.value[0],rangeMaxValue=maxRange.value[1];
		}
		else if (maxRange.type.match(/^qua/)) { // quartile
			//Compute quartiles
			var orderedValues=[];
			for (var i=0; i<rowList.length; i++) {
				for (var j=0; j<rowList[i].heatCellList.length; j++) {
					if (rowList[i].heatCellList[j].value != null) orderedValues.push(rowList[i].heatCellList[j].value);
				}
			}
			orderedValues.sort(sortNumber);
			var lowQuartValue=computeQuartile(orderedValues,orderedValues.length / 4);
			var midQuartValue=computeQuartile(orderedValues,orderedValues.length / 2);
			var upQuartValue=computeQuartile(orderedValues,3 * orderedValues.length / 4);
			rangeMinValue=midQuartValue - (midQuartValue - lowQuartValue) * maxRange.value;
			rangeMaxValue=midQuartValue + (upQuartValue - midQuartValue) * maxRange.value;
//console.log(lowQuartValue,midQuartValue,upQuartValue);
		}
//console.log(rangeMinValue,rangeMaxValue);
		// Update cell values
		for (var i=0; i<rowList.length; i++) {
			for (var j=0; j<rowList[i].heatCellList.length; j++) {
				if (rowList[i].heatCellList[j].value != null) {
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
		if (palette.trimValues.svgMin) {palette.trimValues.svgMin.data('line').remove(); palette.trimValues.svgMin.remove();}
		if (palette.trimValues.svgMax) {palette.trimValues.svgMax.data('line').remove(); palette.trimValues.svgMax.remove();}
		palette.trimValues={min:100,max:100,svgMin:null,svgMax:null};

		if (HM.normProcess.reference=='z-score') {
			if (!maxAbsZscore) {
				processData();
//console.log('Max z-score='+maxAbsZscore);
			}
			/* Convert Z-scores to % */
			//convertZscores2pc();
			for (var i=0; i<rowList.length; i++) {
				for (var j=0; j<rowList[i].heatCellList.length; j++) {
					if (rowList[i].heatCellList[j].value != null && !rowList[i].heatCellList[j].isExcluded) {
						rowList[i].heatCellList[j].pcValue=(rowList[i].heatCellList[j].zscore==null)? 100 : Math.round(100*(rowList[i].heatCellList[j].zscore/maxAbsZscore)); //reference=0
						if (rowList[i].heatCellList[j].pcValue < 0) {palette.hasNegValues=true;}
						else {palette.hasPosValues=true;}
					}
				}
			}
		}
		else {processData();}
		this.computeCellsColor();
//console.timeEnd('COMPUTE');
    }
    function processData() {
	    if (HM.normProcess.scope=='row') {
		    for (var i=0; i<rowList.length; i++) {
			    normalizeData(rowList[i].heatCellList);
		    }
	    }
	    else if (HM.normProcess.scope=='column') {
		    for (var i=0; i<columnList.length; i++) {
			    normalizeData(columnList[i].heatCellList);
		    }
	    }
	    else { // entire data set
		    var cellList=[];
		    for (var i=0; i<rowList.length; i++) {
			    for (var j=0; j<rowList[i].heatCellList.length; j++) {
				    cellList.push(rowList[i].heatCellList[j]);
			    }
		    }
		    normalizeData(cellList);
	    }
    }

    function normalizeData(cellList) {
		var valueList=[];
		for (var i=0; i<cellList.length; i++) {
			if (cellList[i].value != null && !cellList[i].isExcluded) valueList.push(cellList[i].value);
		}
		if (HM.normProcess.reference=='z-score') {
			var mean=0;
			for (var i=0; i<valueList.length; i++) {
				mean+=valueList[i];
			}
			mean/=valueList.length; // mean value
			var sqDelta=0;
			for (var i=0; i<valueList.length; i++) {
				sqDelta+=Math.pow(valueList[i]-mean,2); // delta^2
			}
			var stdDev=Math.sqrt(sqDelta/valueList.length);
//console.log('StdDev='+stdDev);

			for (var i=0; i<cellList.length; i++) { // compute z-scores & record biggest
				if (cellList[i].value != null && !cellList[i].isExcluded) {
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
			valueList.sort(ascNumSort);
			var numValues=valueList.length;
			var minValue=valueList[0];
			var maxValue=valueList[numValues-1];
//for (var i=0; i<valueList.length; i++) {console.log(valueList[i]);}
			var referenceValue=0;
			if (HM.normProcess.reference=='mean') {
				for (var i=0; i<valueList.length; i++) {referenceValue+=valueList[i];}
				referenceValue/=valueList.length;
			}
			if (HM.normProcess.reference=='median') {
				var halfNumValues=numValues/2;
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
			for (var i=0; i<cellList.length; i++) {
				if (cellList[i].value != null) {
					cellList[i].pcValue=Math.round(100*((cellList[i].value-referenceValue)/valueMaxRange));
					if (cellList[i].value < referenceValue) {palette.hasNegValues=true;}
					else {palette.hasPosValues=true;}
				}
			}
//console.log('1',referenceValue,palette.hasNegValues,palette.hasPosValues);
//console.log(HM.normProcess.reference+': Reference='+referenceValue+', Range='+valueMaxRange+' ['+minValue+' -> '+maxValue+'], Num Values='+numValues);

		}
	    //HM.computeCellsColor(cellList);
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
		//cellHighlightColor=(this.normProcess.colors.match('BW'))? '#DD0' : (this.normProcess.colors.match('BY'))? '#0DD' : '#0AA';
		cellHighlightColor=(this.normProcess.colors.match('GB|BB'))? '#FF0' : '#0DD'; //(this.normProcess.colors.match('BW|BY'))? '#0DD' : '#0AA';
		noDataCellColor=(this.normProcess.colors.match('W'))? '#000' : '#EEE';
		if (!cellList) {
//console.log('Cell list...');
			cellList=[];
			for (var i=0; i<rowList.length; i++) {
				for (var j=0; j<rowList[i].heatCellList.length; j++) {
					cellList.push(rowList[i].heatCellList[j]);
				}
			}
		}
//console.log(cellList.length+' cells!');
		var scaleTrimMax=100/palette.trimValues.max;
		var scaleTrimMin=100/palette.trimValues.min;
		for (var i=0; i<cellList.length; i++) {
			var hCell=cellList[i];
			if (hCell.isExcluded) {hCell.color='#999';}
			else if (hCell.value == null) {hCell.color=noDataCellColor;}
			else {
				//if (hCell.value < referenceValue)
				if (hCell.pcValue < 0) { // below reference
					var pcValue=Math.min(100,Math.abs(hCell.pcValue)); //Math.round(100*((referenceValue-hCell.value)/valueMaxRange));
					//var usedValue=(Math.min(100,(100-pcValue)*scaleTrimMin)); // never when trimValues.min=0
					var usedValueUp=(Math.min(100,pcValue*scaleTrimMin));
					var usedValueDown=Math.max(0,100-pcValue*scaleTrimMin);
					//Green
					if (HM.normProcess.colors.match('GB')) {hCell.color='rgb(0%,'+usedValueUp+'%,0%)';} // green->black
					else if (HM.normProcess.colors.match('GW')) {hCell.color='rgb('+usedValueDown+'%,100%,'+usedValueDown+'%)';} // green->white
					if (HM.normProcess.colors.match('GY')) {hCell.color='rgb('+usedValueDown+'%,100%,0%)';} // green->yellow
					//Blue
					if (HM.normProcess.colors.match('BB')) {hCell.color='rgb(0%,0%,'+usedValueUp+'%)';} // blue->black
					else if (HM.normProcess.colors.match('BW')) {hCell.color='rgb('+usedValueDown+'%,'+usedValueDown+'%,100%)';} // blue->white
					if (HM.normProcess.colors.match('BY')) {hCell.color='rgb('+usedValueDown+'%,'+usedValueDown+'%,'+usedValueUp+'%)';} // blue->yellow
				}
				else { // above reference
					var pcValue=Math.min(100,hCell.pcValue); //Math.round(100*((hCell.value-referenceValue)/valueMaxRange));
					var usedValueUp=(Math.min(100,pcValue*scaleTrimMax));
					var usedValueDown=Math.max(0,100-pcValue*scaleTrimMax);
					//Red
					if (HM.normProcess.colors.match('BR')) {hCell.color='rgb('+usedValueUp+'%,0%,0%)';} // black->red
					else if (HM.normProcess.colors.match('WR')) {hCell.color='rgb(100%,'+usedValueDown+'%,'+usedValueDown+'%)';} // white->red
					if (HM.normProcess.colors.match('YR')) {hCell.color='rgb(100%,'+usedValueDown+'%,0%)';} // yellow->red
				}
//console.log(i+': '+hCell.color);
			}
			//hCell.aspect.attr({fill:hCell.color,stroke:hCell.color});
			mapContext.fillStyle=hCell.color;
			//mapContext.globalAlpha = (hCell.isFlagged)? 0.2 : 1;
			mapContext.fillRect(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx,hCellW,hCellH);
			//Flag
			if (hCell.isFlagged && cellFlagVisible) {drawCellFlag(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx);}
			//Selection
			if (hCell.isSelectedData) {hCell.isSelectedData.attr({stroke:cellHighlightColor});}
		}

		drawPalette(noDistrib);

		if (focusArea) {  // update zoom area if already initialized
			focusArea.attr({stroke:cellHighlightColor}); // ,fill:cellHighlightColor
			focusArea.data('close')[0].attr({fill:cellHighlightColor});
			if (focusArea.data('selRows')) {focusArea.data('selRows')[0].attr({fill:cellHighlightColor});}
			if (focusArea.data('selColumns')) {focusArea.data('selColumns')[0].attr({fill:cellHighlightColor});}
		}
	}

	this.updateCellsValue=function(scope,itemId,updateString) {
		var labelIdxList={row:{},column:{}}; // type.id -> idx
		var selLabel;
		for (var type in labelIdxList) {
			var labelList=(type=='row')? rowList : columnList;
			for (var i=0; i<labelList.length; i++) {
				labelIdxList[type][labelList[i].id]=i;
				if (type==scope && labelList[i].id==itemId) {selLabel=labelList[i];}
			}
		}
		var updateData=updateString.split(';');
		for (var i=0; i<updateData.length; i++) {
			var cellData=updateData[i].split(',');
			var rowIdx,colIdx;
			if (scope=='row') {
				rowIdx=labelIdxList.row[itemId];
				colIdx=labelIdxList.column[cellData[0]];
			}
			else { // column
				rowIdx=labelIdxList.row[cellData[0]];
				colIdx=labelIdxList.column[itemId];
			}
			var hCell=rowList[rowIdx].heatCellList[colIdx];
			hCell.value=(cellData[1]==null || cellData[1].length==0 || isNaN(cellData[1]*1))? null : cellData[1]*1; // null*1=''*1=0!!!!
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
	}

	function updateCellsExclusion(label,excludedStatus) { // label object,true/false
		for (var i=0; i<label.heatCellList.length; i++) {
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
			var binSize,halfRange;
			if (palette.hasNegValues && palette.hasPosValues) {
				binSize=5;
				halfRange=100/binSize;
			}
			else {
				binSize=2.5;
				halfRange=(palette.hasNegValues)? 100/binSize : 0;
			}
			var binHalf=binSize/2;
			for (var b=0; b<=40; b++) {palette.valueDistrib[b]=0;} // 41 bins
			for (var i=0; i<rowList.length; i++) {
				for (var j=0; j<rowList[i].heatCellList.length; j++) {
					var hCell=rowList[i].heatCellList[j];
					if (hCell.isExcluded || hCell.value==null) {continue;}
					var bin;
					if (hCell.pcValue < 0) {
						var pcValue=Math.min(100,Math.abs(hCell.pcValue));
						bin=halfRange-Math.floor((pcValue+binHalf)/binSize);
					}
					else {
						var pcValue=Math.min(100,hCell.pcValue);
						bin=halfRange+Math.floor((pcValue+binHalf)/binSize);
					}
					palette.valueDistrib[bin]++;
				}
			}
		}
//console.log(palette.valueDistrib);

		/* (Re)Draw palette */
		var palSVG=palette.svgElements,palDistr=palette.valueDistrib;
		var palX=palette.posX+0.5, palY=palette.posY+0.5;
		//Clear palette
		if (palSVG.length) {
			for (var i=0; i<palSVG.length; i++) {palSVG[i].remove();}
		}
		palSVG.length=0;

		var palH=100,palW=200; // palette dimensions
		//Draw color gradients
		var colorsL=HM.normProcess.colorList.limit;
		var colorsR=HM.normProcess.colorList.reference;
		var grad=HM.normProcess.colors.split('');
		if (palette.hasNegValues && palette.hasPosValues) {
			//Negative values
			palSVG.push(HM.paper.rect(palX+100-palette.trimValues.min,palY,palette.trimValues.min,100).attr({stroke:'none',fill:'180-'+colorsR[grad[1]]+'-'+colorsL[grad[0]]})); // gradient
			palSVG.push(HM.paper.rect(palX,palY,100-palette.trimValues.min,100).attr({stroke:'none',fill:colorsL[grad[0]]}));
			//Positive values
			palSVG.push(HM.paper.rect(palX+100,palY,palette.trimValues.max,100).attr({stroke:'none',fill:'0-'+colorsR[grad[1]]+'-'+colorsL[grad[2]]})); // gradient
			palSVG.push(HM.paper.rect(palX+99+palette.trimValues.max,palY,100-palette.trimValues.max,100).attr({stroke:'none',fill:colorsL[grad[2]]}));
			palette.refLimit=palX+100;
		}
		else {
			if (palette.hasNegValues) {
				var range=palette.trimValues.min*2;
				palSVG.push(HM.paper.rect(palX+200-range,palY,range,100).attr({stroke:'none',fill:'180-'+colorsR[grad[1]]+'-'+colorsL[grad[0]]}));
				palSVG.push(HM.paper.rect(palX,palY,200-range,100).attr({stroke:'none',fill:colorsL[grad[0]]}));
				palette.refLimit=palX+200;
			}
			else {
				var range=palette.trimValues.max*2;
				palSVG.push(HM.paper.rect(palX,palY,range,100).attr({stroke:'none',fill:'0-'+colorsR[grad[1]]+'-'+colorsL[grad[2]]}));
				palSVG.push(HM.paper.rect(palX+range,palY,200-range,100).attr({stroke:'none',fill:colorsL[grad[2]]}));
				palette.refLimit=palX;
			}
		}
//console.log(palette.trimValues.min,palette.trimValues.max);
		//Daw histogram of distribution
		var binMaxValue=0;
		for (var i=0; i<palDistr.length; i++) {binMaxValue=Math.max(binMaxValue,palDistr[i]);}
		var palScale=palH/(1.1*binMaxValue);
		var dpStrg='M'+palX+','+(palY+palH);
		var curY=palY;
		for (var i=0; i<=40; i++) { // 41 bins
			var dx=(i==0 || i==40)? 3 : 5;
			var y=palY-Math.round(palDistr[i]*palScale);
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
		var trim=palette.trimValues;
		if (palette.hasNegValues) {
			if (trim.svgMin) {trim.svgMin.data('line').attr({stroke:cellHighlightColor});}
			else {
				var tX=palX+(100-trim.min)*palScale;
				var line=HM.paper.path('M'+tX+','+palY+' l0,100').attr({stroke:cellHighlightColor,'stroke-width':2}).hide();
				trim.svgMin=HM.paper.path('M'+tX+','+(palY+101)+' l8,8 l-17,0 Z').attr({fill:'#000',stroke:'none'})
				.data({type:'min',pos:tX,totShift:0,shift:0,line:line}).hover(function(){this.attr({fill:'#F00'});},function(){if(!dragContext){this.attr({fill:'#000'})}})
				.drag(dragTrimMove,dragTrimStart,dragTrimStop);
			}
		}
		if (palette.hasPosValues) {
			if (trim.svgMax) {trim.svgMax.data('line').attr({stroke:cellHighlightColor});}
			else {
				var tX=palX+palW-(100-trim.max)*palScale;
				var line=HM.paper.path('M'+tX+','+palY+' l0,100').attr({stroke:cellHighlightColor,'stroke-width':2}).hide();
				trim.svgMax=HM.paper.path('M'+tX+','+(palY+101)+' l8,8 l-17,0 Z').attr({fill:'#000',stroke:'none'})
				.data({type:'max',pos:tX,totShift:0,shift:0,line:line}).hover(function(){this.attr({fill:'#F00'});},function(){if(!dragContext){this.attr({fill:'#000'})}})
				.drag(dragTrimMove,dragTrimStart,dragTrimStop);
			}
		}
    }


    function ascNumSort(a,b) {return a-b}


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
		var mousePos=getMousePositionInElement(e); // relative to hCellsSensor NOT chart (not clear why; due to multiple DIVs?)

		var colIdx=Math.round((mousePos[0]/hCellW)-0.5); // - because index not rank
		var rowIdx=Math.round((mousePos[1]/hCellH)-0.5);
//document.getElementById('mX').value=rowIdx;
//document.getElementById('mY').value=colIdx;
		if (rowIdx>=rowList.length) {rowIdx=rowList.length-1;}
		else if (rowIdx<0) {rowIdx=0;}
		if (colIdx>=columnList.length) {colIdx=columnList.length-1;}
		else if (colIdx<0) {colIdx=0;}
		return rowList[rowIdx].heatCellList[colIdx];
	}
    function setCellEmphasis(hCell,action,e) {
		var cellBorderColor,cellCursor; //labelColor,
		//var c=hCell.aspect;
		var cX=hCellW*hCell.columnIdx, cY=hCellH*hCell.rowIdx;
		var rowLabel=rowList[hCell.rowIdx],colLabel=columnList[hCell.columnIdx];
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
			var valueStrg;
			if (cellOnMouseOver) {valueStrg=cellOnMouseOver(hCell.value,hCell,rangeMinValue,rangeMaxValue);} // hCell in case more complex display is needed
			else if (hCell.value==null) {valueStrg='* no value *';}
			else {
				var overStrg=(hCell.value <= rangeMinValue)? '≤' : (hCell.value >= rangeMaxValue)? '≥' : '';
				valueStrg='Value: ' + overStrg + (Math.round(100*hCell.value)/100);
				if (HM.normProcess.reference=='z-score' && !hCell.isExcluded) {
					if (hCell.zscore==null) {valueStrg+=' [no z-score]';}
					else {valueStrg+=' [z='+(Math.round(100*hCell.zscore)/100)+']';}
				}
			}
			var textStrg=(e.ctrlKey || e.shiftKey)? HM.entities.row+': '+rowLabel.label+'\n'+HM.entities.column+': '+colLabel.label+'\n'+valueStrg : valueStrg;
			if (hCell.isFlagged) textStrg+='\n* '+flagText+' *';
			//var textStrg=HM.entities.row+': '+rowLabel.label+'\n'+HM.entities.column+': '+colLabel.label+'\n'+valueStrg;
//textStrg+='\n'+hCell.rowIdx+','+hCell.columnIdx;
			var x=cellX0+cX+Math.round(hCellW/2); // c.attr('x')+Math.round(c.attr('width')/2);
			var y=cellY0+cY; //c.attr('y');
			drawCellPopup(x,y,textStrg,cellHighlightColor,cellHighlightBgColor,cellHighlightBgColor);
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
		}
		//c.attr({stroke:cellBorderColor,cursor:cellCursor});
		mapContext.strokeStyle=cellBorderColor;
		mapContext.lineWidth=2;
		mapContext.strokeRect(cX+1,cY+1,hCellW-2,hCellH-2);
		if (hCell.isFlagged && action=='off' && cellFlagVisible) {drawCellFlag(cX,cY);}
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
		var t=HM.paper.text(x,y,textStrg).attr({'font-size':10,'text-anchor':'start',fill:textColor});
		var tBox0=t.getBBox();
		t.attr({x:x-Math.round(tBox0.width/2),y:y-Math.round(tBox0.height/2)-5});
		var tBox=t.getBBox();
		var tx=tBox.x;
		var ty=tBox.y;
		var tw=tBox.width;
		var th=tBox.height;
		var tb=HM.paper.rect(tx-2,ty-1,tw+4,th+2,4).attr({stroke:borderColor,fill:backColor,'fill-opacity':0.75});
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
			var selStatus;
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
			var midX=cellX0 + (hCell.columnIdx * hCellW) + Math.round(hCellW / 2);
			var bottomY=cellY0 + (hCell.rowIdx * hCellH) + hCellH;
			cellOnClick(rowList[hCell.rowIdx].id,columnList[hCell.columnIdx].id,selStatus,midX,bottomY); // user callback (rowId,colId,selStatus,midX,bottomY)
		}
	}
	function drawCellSelection(hCell,shiftX,shiftY) {
		if (!shiftX) shiftX=0;
		if (!shiftY) shiftY=0;
		var s;
		var cX=cellX0 + (hCell.columnIdx * hCellW) + Math.round(hCellW / 2) + shiftX;
		var cY=cellY0 + (hCell.rowIdx * hCellH) + Math.round(hCellH / 2) + shiftY;
		if (hCellW >=20 && hCellH >=20) {
			var x= cX - 3;
			var y= cY - 2;
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
			var label=(HM.normProcess.scope=='row')? rowList[hCell.rowIdx] : columnList[hCell.columnIdx];
			normalizeData(label.heatCellList);
			HM.computeCellsColor(label.heatCellList);
		}
		//HM.colorCells();
		if (cellOnModClick) {
			var midX=cellX0 + (hCell.columnIdx * hCellW) + Math.round(hCellW / 2);
			var bottomY=cellY0 + (hCell.rowIdx * hCellH) + hCellH;
			cellOnModClick(rowList[hCell.rowIdx].id,columnList[hCell.columnIdx].id,hCell.isExcluded,midX,bottomY); // user callback (rowId,colId,selStatus,midX,bottomY)
		}
    }
	this.clearAllCells=function() {
		for (var i=0; i<rowList.length; i++) {
			for (var j=0; j<rowList[i].heatCellList.length; j++) {
				var hCell=rowList[i].heatCellList[j];
				if (hCell.isSelectedData) {
					hCell.isSelectedData.remove();
					hCell.isSelectedData=null;
				}
			}
		}
		if (oneCellSelection) {lastSelectedCell=null;}
	}
	this.flagCell=function(colLabelID,rowLabelID,status) {
		//Finding cell coordinates
		var colIdx=-1, rowIdx=-1;
		for (var i=0; i<columnList.length; i++) {
			if (columnList[i].id==colLabelID) {
				colIdx=i;
				break;
			}
		}
		if (colIdx==-1) {return -1;}
		for (var i=0; i<rowList.length; i++) {
			if (rowList[i].id==rowLabelID) {
				rowIdx=i;
				break;
			}
		}
		if (rowIdx==-1) {return -2;}
		var hCell=rowList[rowIdx].heatCellList[colIdx];
		hCell.isFlagged=status;
		var color=(status==true)? cellHighlightColor : (hCell.isExcluded)? '#999' : hCell.color;
		if (cellFlagVisible) drawCellFlag(hCellW*colIdx,hCellH*rowIdx,color);
		return 0;
	}
	this.showCellFlag=function(visStatus) {
		cellFlagVisible=(visStatus==false)? false : true;
		for (var i=0; i<rowList.length; i++) {
			for (var j=0; j<rowList[i].heatCellList.length; j++) {
				if (rowList[i].heatCellList[j].isFlagged) {
					var hCell=rowList[i].heatCellList[j];
					//canvas
					if (visStatus==false) { // redraw cell w/o flag
						mapContext.fillStyle=hCell.color;
						mapContext.fillRect(hCellW*hCell.columnIdx,hCellH*hCell.rowIdx,hCellW,hCellH);
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
	}

	/**********************************************************/
    /********************* Palette events *********************/
    /**********************************************************/
	function dragTrimStart() {
		dragContext=true;
		this.data('line').toFront().show();
	}
	function dragTrimMove(dx,dy) {
		var midLimit=palette.refLimit;
		var shift=this.data('totShift')+dx;
		var tX=this.data('pos')+shift;
		if (this.data('type')=='min') {
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
		var scale=(palette.hasNegValues && palette.hasPosValues)? 1 : 2;
		if (this.data('type')=='min') {
			palette.trimValues.min=Math.round((palette.refLimit-this.data('pos')-this.data('totShift'))/scale);
		}
		else {
			palette.trimValues.max=Math.round((this.data('pos')+this.data('totShift')-palette.refLimit)/scale);
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
		var type=(l.data('row'))? 'row' : 'column';
		var returnAction;
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
			var excludedStatus=(returnAction=='exclude')? true : false;
			updateCellsExclusion(l.data(type),excludedStatus);
		}
	}
	function labelDoubleClick(l,e) {
		dblClickContext=true;
//console.log('*DBL-CLICK '+dblClickContext);
		var type=(l.data('row'))? 'row' :'column';
		var itemLabel=l.data(type);
		var lBox=l.getBBox();
		if (type=='row') {
			var posX=(rowLabelLocation=='left')? lBox.x2 : lBox.x;
			if (e.ctrlKey || e.shiftKey) {if (rowOnModDblClick) {rowOnModDblClick(itemLabel.id,posX,lBox.y);}}
			else if (rowOnDblClick) {rowOnDblClick(itemLabel.id,posX,lBox.y);}
		}
		else { // column
			var posY=(colLabelLocation=='top')? lBox.y2 : lBox.y;
			if (e.ctrlKey || e.shiftKey) {if (columnOnModDblClick) {columnOnModDblClick(itemLabel.id,lBox.x2,posY);}}
			else if (columnOnDblClick) {columnOnDblClick(itemLabel.id,lBox.x2,posY);}
		}
		// Select all cells in row or column
		if (cellOnClick && !oneCellSelection) {
			var newSelStatus=(e.ctrlKey || e.shiftKey)? false : true;
			for (var i=0; i<itemLabel.heatCellList.length; i++) {
				var hCell=itemLabel.heatCellList[i];
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
		var itemLabel,type,isLabel=false,isAnnot=false,isGroup=false,typePopup;
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
		var labelColor,labelCursor;
		if (action=='on') {
			if (labelPopupText) return;
			labelColor='#F00';
			labelCursor=(!isLabel)? 'pointer' : (movableLabel[type])? 'move' : 'auto';
			var popupText=itemLabel.popupText || itemLabel.label
			//if (itemLabel.popupText) {
				if (typePopup=='row') {
					var x=(rowLabelLocation=='left')? l.attr('x')+6 : l.attr('x')-6; //mapGeometry[rowLabelLocation];
					var y=l.attr('y');
					//drawPopup(x,y,itemLabel.popup,'#000','#FFF');
					var t=l.paper.text(x,y,popupText).attr({'font-size':12,fill:'#000','text-anchor':'start'});
					var tBox=t.getBBox();
					if (rowLabelLocation=='right') {t.attr('x',t.attr('x')-tBox.width);}
					if (tBox.y < 0) {t.attr('y',t.attr('y')-tBox.y+4);} // tBox.y is negative
					else if (tBox.y+tBox.height > paperHeight) {t.attr('y',t.attr('y')-((tBox.y+tBox.height)-paperHeight+4));}
					tBox=t.getBBox();
					var tx=tBox.x;
					var ty=tBox.y;
					var tw=tBox.width;
					var th=tBox.height;
					var tb=l.paper.rect(tx-2,ty-1,tw+4,th+2,4).attr({stroke:'#000',fill:'#FFF','fill-opacity':1});
					t.toFront();
					labelPopupText=l.paper.set(tb,t);
					labelPopupText.transform(l.transform()); // in case inside focusArea
				}
				else { // column
					var x=l.attr('x');
					var y=(colLabelLocation=='top')? l.attr('y')+12 : l.attr('y')-12; //mapGeometry[rowLabelLocation];
					var t=l.paper.text(x,y,popupText).attr({'font-size':12,fill:'#000','text-anchor':'middle'});
					var tBox=t.getBBox();
					if (tBox.x < 0) {t.attr('x',t.attr('x')-tBox.x+4);} // tBox.x is negative
					else if (tBox.x+tBox.width > paperWidth) {t.attr('x',t.attr('x')-((tBox.x+tBox.width)-paperWidth+4));}
					tBox=t.getBBox();
					var tx=tBox.x;
					var ty=tBox.y;
					var tw=tBox.width;
					var th=tBox.height;
					var tb=l.paper.rect(tx-2,ty-1,tw+4,th+2,4).attr({stroke:'#000',fill:'#FFF','fill-opacity':1});
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
				var labelList=(type=='row')? rowList : columnList;
				for (var i=itemLabel.labelIndexRange[0]; i<=itemLabel.labelIndexRange[1]; i++) {
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
				var labelList=(type=='row')? rowList : columnList;
				for (var i=itemLabel.labelIndexRange[0]; i<=itemLabel.labelIndexRange[1]; i++) {
					labelList[i].aspect.attr({fill:colorList[labelList[i].selectColorIdx]});
				}
			}
		}
		l.attr({fill:labelColor,cursor:labelCursor});
    }

	/* Multi-labels external action */
	function focussedLabelsAction(call,type) {
		var labelList,selLabelList;
		if (call=='tree') {
			labelList=(type=='row')? columnList : rowList; //row/col inverted for tree
			selLabelList=selectedLabels;
			if (labelsOnSelect[type].ranking) {
				var rankedData=[];
				for (var i=0; i<labelList.length; i++) {rankedData.push(labelList[i].id);}
				labelsOnSelect[type].action(selLabelList,rankedData);
			}
			else {labelsOnSelect[type].action(selLabelList);}
		}
		else { // call=focus
			var targetList,selIdx,idx;
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
			for (var i=focusArea.data('labelIdx')[selIdx]; i<=focusArea.data('labelIdx')[selIdx+1]; i++) {selLabelList.push(targetList[i].id);}
			if (labelsOnSelect[type].ranking) {
				var rankedData=[];
				for (var i=focusArea.data('labelIdx')[idx]; i<=focusArea.data('labelIdx')[idx+1]; i++) {rankedData.push(labelList[i].id);}
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
			for (var i=0; i<rowList.length; i++) {
				dragParams.labelList[rowList[i].rank]=rowList[i];
			}
			dragParams.maxRank=rowList.length;
			dragParams.type='row';
			dragParams.axis='y';
			dragParams.emptyPosL=label.attr('y');
		}
		else {
			dragParams.itemLabel=label.data('column'); // row or column label
			for (var i=0; i<columnList.length; i++) {
				dragParams.labelList[columnList[i].rank]=columnList[i];
			}
			dragParams.maxRank=columnList.length;
			dragParams.type='column';
			dragParams.axis='x';
			dragParams.emptyPosL=label.attr('x');
		}
		dragParams.startPosL=dragParams.emptyPosL; // original pos at drag start
		dragParams.dragPrevLength=0;
		label.toFront();
		label.attr('fill','#F00');
		dragParams.svgCells=[];

		if (dragParams.itemLabel.heatCellList.length <= 2) { // few cells to move
			for (var i=0; i<dragParams.itemLabel.heatCellList.length; i++) {
				var hCell=dragParams.itemLabel.heatCellList[i];
				var svgCell=HM.paper.rect(cellX0+hCellW*hCell.columnIdx+1,cellY0+hCellH*hCell.rowIdx+1,hCellW-2,hCellH-2,0).attr({fill:hCell.color,'fill-opacity':0.5,'stroke-width':2,stroke:cellHighlightColor});
				dragParams.svgCells.push(svgCell);
				if (hCell.isSelectedData) {hCell.isSelectedData.toFront();}
			}
		}
		else { // draw a single svg rect for all cells
			var setX,setY,setW,setH; // row/column of cells selected for drag
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
			var svgCellSet=HM.paper.rect(setX,setY,setW,setH,0).attr({fill:cellHighlightColor,'fill-opacity':0.5,stroke:'none'});
			dragParams.svgCells.push(svgCellSet);
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
			for (var i=0; i<dragParams.svgCells.length; i++) {
				dragParams.svgCells[i].remove();
			}
			//setTimeout(function(){dragContext=null; dragParams={};},300); // delay set dragContext to null to prevent click action
		}
		dragParams={};
	}
    function dragLabelMove(dx,dy) {
		if (!dragContext) { // 1st move -> initialize dragParams
			dragLabelInit(this);
		}
		if (labelPopupText) {labelPopupText.remove(); labelPopupText=null;}
		dragContext=true;
		var dragLength;
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
			var tr=this.transform().toString().replace(/r[^t]+/,''); // in case inside focusArea
			//this.transform(''); // unrotate
			this.attr(dragParams.axis,dragParams.startPosL+dragLength);
			//this.rotate(colLabelAngle,this.attr('x'),this.attr('y')); // rerotate
			this.transform(tr+'r'+colLabelAngle+','+this.attr('x')+','+this.attr('y')); // rerotate
		}
		/* Moving svg cells */
		var dragDelta=dragLength-dragParams.dragPrevLength;
		for (var i=0; i<dragParams.svgCells.length; i++) {
			dragParams.svgCells[i].attr(dragParams.axis,dragParams.svgCells[i].attr(dragParams.axis)+dragDelta);
		}
		/* Moving cell selection (if any) */
		if (dragParams.svgCells.length > 1) { // tempory svg cells
			var shiftX=0,shiftY=0;
			if (dragParams.axis=='x') {shiftX=this.attr(dragParams.axis)-dragParams.emptyPosL;}
			else {shiftY=this.attr(dragParams.axis)-dragParams.emptyPosL;}
			for (var i=0; i<this.data(dragParams.type).heatCellList.length; i++) {
				var hCell=this.data(dragParams.type).heatCellList[i];
				if (hCell.isSelectedData) {
					hCell.isSelectedData.remove();
					drawCellSelection(hCell,shiftX,shiftY);
				}
			}
		}
		/* moving annotation cells */
		//var shift=dragParams.emptyPosL-this.attr(dragParams.axis);
		for (var setName in dragParams.itemLabel.annotations) {
		    var aC=dragParams.itemLabel.annotations[setName].aspect;
			aC.attr(dragParams.axis,aC.attr(dragParams.axis)+dragDelta);
		}
		/* Checking if dragging over another label */
		if (dragDelta < 0) { // dragging up/left
			if (dragParams.itemLabel.rank > 1) {
				var overLabel=dragParams.labelList[dragParams.itemLabel.rank-1];
				if (this.attr(dragParams.axis)-overLabel.aspect.attr(dragParams.axis) < maxOverlap[dragParams.type]) {
					clearCanvasRowColumn(overLabel,dragParams.type);
					moveLabelAndCells(overLabel.aspect); // false => move canvas cells
				}
			}
		}
		else { // dragging down/right
			if (dragParams.itemLabel.rank < dragParams.maxRank) { // max rank
				var overLabel=dragParams.labelList[dragParams.itemLabel.rank+1];
				if (overLabel.aspect.attr(dragParams.axis)-this.attr(dragParams.axis) < maxOverlap[dragParams.type]) { //-2 for safety
					clearCanvasRowColumn(overLabel,dragParams.type);
					moveLabelAndCells(overLabel.aspect); // false => move canvas cells
				}
			}
		}
		dragParams.dragPrevLength=dragLength;
    }
	function clearCanvasRowColumn(label,type) {
		var setX,setY,setW,setH; // row/column of cells selected for drag
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
		var prevPosL=movedLabel.attr(dragParams.axis);
		/* moving canvas cells */
		var shift=dragParams.emptyPosL-prevPosL;
		moveCells(movedLabel.data(dragParams.type),dragParams.axis,shift,endDrag);
		/* moving label */
		var dragList,otherList,dragTypeIdx;
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
		var overRank=movedLabel.data(dragParams.type).rank;
		var itemRank=dragParams.itemLabel.rank;
		if (overRank != itemRank) { // not end of dragging
			movedLabel.data(dragParams.type).rank=itemRank;
			dragParams.itemLabel.rank=overRank;
			dragParams.labelList[itemRank]=movedLabel.data(dragParams.type);
			dragParams.labelList[overRank]=dragParams.itemLabel;
			//indexes
			var overCellList=dragList[overRank-1];
			dragList[overRank-1]=dragList[itemRank-1];
			dragList[itemRank-1]=overCellList;
			for (var i=0; i<otherList.length; i++) {
				var overCell=otherList[i].heatCellList[overRank-1];
				otherList[i].heatCellList[overRank-1]=otherList[i].heatCellList[itemRank-1];
				otherList[i].heatCellList[itemRank-1]=overCell;
				otherList[i].heatCellList[overRank-1][dragTypeIdx]=overRank-1;
				otherList[i].heatCellList[itemRank-1][dragTypeIdx]=itemRank-1;
			}
		}
		/* Label rotation & translation if inside focusArea */
		var overIdx=itemRank-1; // label were switched or end of drag
		if (dragParams.type=='row') {
			var trF='';
			if (focusArea.data('status')=='on' && overIdx >= focusArea.data('labelIdx')[0] && overIdx <= focusArea.data('labelIdx')[1]) { // inside focus zone
				var focusShift=(rowLabelLocation=='left')? focusArea.attr('x')-cellX0 : focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
				trF='t'+focusShift+',0';
				movedLabel.toFront(); // shift & move above bgRow
			}
			movedLabel.transform(trF); // also when empty in case moved outside focus
		}
		else {
			var trF='';
			if (focusArea.data('status')=='on' && overIdx >= focusArea.data('labelIdx')[2] && overIdx <= focusArea.data('labelIdx')[3]) { // inside focus zone
				var focusShift=(colLabelLocation=='top')? focusArea.attr('y')-cellY0 : focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
				trF='t0,'+focusShift;
				movedLabel.toFront(); // move above bgColumn
			}
			movedLabel.transform(trF+'r'+colLabelAngle+','+movedLabel.attr('x')+','+movedLabel.attr('y'));
		}

		/* storing new empty pos */
		dragParams.emptyPosL=prevPosL;
    }

    function moveCells(label,axis,shift,endDrag) { // axis: x (column) or y (row)!
		var cList=label.heatCellList;
		var shiftX,shiftY;
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
		for (var i=0; i < cList.length; i++) {
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
		for (var setName in label.annotations) {
		    var aC=label.annotations[setName].aspect;
			aC.attr(axis,aC.attr(axis)+shift);
		}
    }

	this.renameLabel=function(type,labelID,newName) {
		var labelList=(type=='row')? rowList : columnList;
		for (var i=0; i<labelList.length; i++) {
			if (labelList[i].id==labelID) {
				labelList[i].label=newName;
				labelList[i].aspect.attr('text',newName);
				break;
			}
		}
	}


	/**********************************************************/
    /************** Drag-drawing zoom area events *************/
    /**********************************************************/
	function dragFocusStart(x,y,e) { // x,y are absolute to document // *public*
		if (focusArea.data('status') == 'on') return; // because triggered by simple click (prevents reset of focusArea on simple click)
		//var mouseInChart=getMousePositionInChart(e);
//debug('ctrl='+modKeyPressed+' : '+x+'('+mousePos[0]+'), '+y+'('+mousePos[1]+')');
		var mousePos=getMousePositionInElement(e);
		var mouseX=mousePos[0]+cellX0,mouseY=mousePos[1]+cellY0;
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
dragContext=true; // ?
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
//console.log('DRAG_1');
		if (dx > 0) {
			focusArea.attr({width:dx});
			if (focusArea.attr('x')+focusArea.attr('width') > mapW+cellX0) {focusArea.attr({width:mapW-focusArea.attr('x')+cellX0})};
		}
		else {
			focusArea.attr({x:focusArea.data('startX')+dx,width:-dx});
			if (focusArea.attr('x') < cellX0) {
				var extra=cellX0 - focusArea.attr('x');
				focusArea.attr({x:cellX0,width:focusArea.attr('width')-extra});
			}
		}
		if (dy > 0) {
			focusArea.attr({height:dy});
			if (focusArea.attr('y')+focusArea.attr('height') > mapH+cellY0) {focusArea.attr({height:mapH-focusArea.attr('y')+cellY0})};
		}
		else {
			focusArea.attr({y:focusArea.data('startY')+dy,height:-dy});
			if (focusArea.attr('y') < cellY0) {
				var extra=cellY0 - focusArea.attr('y');
				focusArea.attr({y:cellY0,height:focusArea.attr('height')-extra});
			}
		}
	}
	function dragFocusEnd() { // *public*
		if (focusArea.data('status') == 'on') {return;} // because triggered by simple click
//console.log('DRAG_2');
		dragContext=false;
		if (focusArea.attr('width') >= hCellW && focusArea.attr('height') >= hCellH) { // big enough => end extension
			focusArea.data({status:'on'});
			/* Find matching row & column label ranges */
			var colIdx1=Math.max(0,Math.round(((focusArea.attr('x')-cellX0)/hCellW)-0.5));
			var colIdx2=Math.min(columnList.length-1,Math.round(((focusArea.attr('x')+focusArea.attr('width')-cellX0)/hCellW)-0.5));
			var rowIdx1=Math.max(0,Math.round(((focusArea.attr('y')-cellY0)/hCellH)-0.5));
			var rowIdx2=Math.min(rowList.length-1,Math.round(((focusArea.attr('y')+focusArea.attr('height')-cellY0)/hCellH)-0.5));
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
			var numRowAnnot=0,numColAnnot=0;
			for (var sName in annotationList.row) {numRowAnnot++;}
			for (var sName in annotationList.column) {numColAnnot++;}
			var annotSpaceX=numRowAnnot*annotationCellSize.row;
			var annotSpaceY=numColAnnot*annotationCellSize.column;
			var cX=(rowLabelLocation=='left')? focusArea.attr('x')+focusArea.attr('width')+annotSpaceX-6 : focusArea.attr('x')-annotSpaceX-28;
			var setY=focusArea.attr('y')-6;
			focusArea.data('close').transform('T'+cX+','+setY).show().toFront();
			if (focusArea.data('selRows')) {focusArea.data('selRows').transform('T'+cX+','+setY).show().toFront();}
			if (focusArea.data('selColumns')) {focusArea.data('selColumns').transform('T'+cX+','+setY).show().toFront();}

			/* Move labels & annotations around focus area */
			//Rows
			var bgW=0;
			for (var i=rowIdx1; i<=rowIdx2; i++) {bgW=Math.max(bgW,rowList[i].aspect.getBBox().width);}
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
			var shiftX,shiftX2;
			if (rowLabelLocation=='left') {
				shiftX=focusArea.attr('x')-cellX0;
				shiftX2=focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
			}
			else {
				shiftX=focusArea.attr('x')+focusArea.attr('width')-cellX0-mapW;
				shiftX2=focusArea.attr('x')-cellX0;
			}
			for (var i=rowIdx1; i<=rowIdx2; i++) {
				rowList[i].aspect.toFront().transform('t'+shiftX+',0');
				/* Moving annotation cells */
				for (var setName in rowList[i].annotations) {
					rowList[i].annotations[setName].aspect.toFront().transform('t'+shiftX2+',0');
				}
			}
			//Columns
			var bgH=0;
			for (var i=colIdx1; i<=colIdx2; i++) {bgH=Math.max(bgH,columnList[i].aspect.getBBox().height);}
			bgH+=10;
			if (colLabelLocation=='top') {
				var shiftD=(colLabelOrientation=='diagonal')? columnList[colIdx2].aspect.getBBox().width : 0;
				focusArea.data('bgColumn').attr({x:focusArea.attr('x')-2,y:focusArea.attr('y')-bgH,width:focusArea.attr('width')+shiftD+4,height:bgH});
				focusArea.data('bgColumn2').attr({x:focusArea.attr('x'),y:focusArea.attr('y')+focusArea.attr('height')+2,width:focusArea.attr('width'),height:annotSpaceY});
			}
			else {
				var shiftD=(colLabelOrientation=='diagonal')? columnList[colIdx1].aspect.getBBox().width : 0;
				focusArea.data('bgColumn').attr({x:focusArea.attr('x')-shiftD-2,y:focusArea.attr('y')+focusArea.attr('height'),width:focusArea.attr('width')+shiftD+4,height:bgH});
				focusArea.data('bgColumn2').attr({x:focusArea.attr('x'),y:focusArea.attr('y')-annotSpaceY-2,width:focusArea.attr('width'),height:annotSpaceY});
			}
			focusArea.data('bgColumn').show().toFront();
			focusArea.data('bgColumn2').show().toFront();
			var shiftY,shiftY2;
			if (colLabelLocation=='top') {
				shiftY=focusArea.attr('y')-cellY0;
				shiftY2=focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
			}
			else {
				shiftY=focusArea.attr('y')+focusArea.attr('height')-cellY0-mapH;
				shiftY2=focusArea.attr('y')-cellY0;
			}
			for (var i=colIdx1; i<=colIdx2; i++) {
				var tr0=columnList[i].aspect.transform();
				columnList[i].aspect.toFront().transform('t0,'+shiftY+tr0);
//console.log('DRAG: '+columnList[i].aspect.transform());  //.toString()
				/* moving annotation cells */
				for (var setName in columnList[i].annotations) {
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
		for (var i=focusArea.data('labelIdx')[0]; i<=focusArea.data('labelIdx')[1]; i++) {
			rowList[i].aspect.transform('');
			for (var setName in rowList[i].annotations) { // annotation cells
				rowList[i].annotations[setName].aspect.transform('');
			}
		}
		//Columns: remove translation, keep original rotation
		for (var i=focusArea.data('labelIdx')[2]; i<=focusArea.data('labelIdx')[3]; i++) {
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
		if (type=='row' && rowList.length==0) {return true}; // no rows => 1D-vertical clustering BUT 'true'=silent validation
		//Check if exists already
		if (annotationList[type][setName]) {
			alert('Annotation \''+setName+'\' already displayed!');
			return false;
		}
		var startShift=annotationCellSize[type]+2;
		var minSpace=(branchList[type])? 50 : -50;
		if (treeSpace[type]-startShift < minSpace) {
			alert('Sorry: No space left. Please delete 1 annotation first.');
			return false;
		}
		closeFocusArea(); // to prevent dealing with its update
		treeSpace[type]-=startShift;
		var labelList;
		var annotStart=treeStart[type];
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

	    var labelIdList={}; // id -> Label (list of available ids)
	    for (var i=0; i<labelList.length; i++) {
			labelIdList[labelList[i].id]=labelList[i];
	    }
		var annotationSet=new AnnotationSet(type,setName);
		var l; // set label SVG object
		if (type=='row') { // Row annot name written alongside the column labels
			var lX=annotStart+Math.round(annotationCellSize[type]/2);
			var lY,lAnchor;
			//if (colLabelLocation=='top') {
				lY=cellY0-5; // always top not to interfere with legends (PP 02/10/16)
				lAnchor='start';
			//}
			//else {
			//	lY=cellY0+mapH+5;
			//	lAnchor='end';
			//}
			//var dispName=(setName.length <= maxColLabelSize)? setName : shortenText(setName,maxColLabelSize);
			//l=HM.paper.text(lX,lY,dispName).attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':lAnchor}).transform('r'+colLabelAngle+','+lX+','+lY);
			l=HM.paper.text(lX,lY,setName).attr({'font-size':colLabelFontSize,'font-weight':'bold','text-anchor':lAnchor});
			var fontSize=colLabelFontSize;
			while (l.getBBox().width*colLabelRatio > colLabelSpace) {
				if (fontSize <= 6) {
					var ratio=colLabelSpace/l.getBBox().width*colLabelRatio;
					var dispName=shortenText(setName,Math.round(setName.length*ratio));
					l.attr({'text':dispName});
					break;
				};
				fontSize-=2;
				l.attr({'font-size':fontSize});
			}
			var angle=(rowLabelLocation=='left' && colLabelLocation=='top')? colLabelAngle : -90; // -90 no diagonal to prevent colision with tree/annot cell
			l.transform('r'+angle+','+lX+','+lY);
		}
		else { // Column annot name written alongside the row labels
			var lY=annotStart+Math.round(annotationCellSize[type]/2);
			var lX,lAnchor;
			if (rowLabelLocation=='left') {
				lX=cellX0-5;
				lAnchor='end';
			}
			else {
				lX=cellX0+mapW+5;
				lAnchor='start';
			}
			//var dispName=(setName.length <= maxRowLabelSize)? setName : shortenText(setName,maxRowLabelSize);
			//l=HM.paper.text(lX,lY,dispName).attr({'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':lAnchor});
			l=HM.paper.text(lX,lY,setName).attr({'font-size':rowLabelFontSize,'font-weight':'bold','text-anchor':lAnchor});
			var fontSize=rowLabelFontSize;
			while (l.getBBox().width > rowLabelSpace) {
				if (fontSize <= 6) {
					var ratio=rowLabelSpace/l.getBBox().width;
					var dispName=shortenText(setName,Math.round(setName.length*ratio));
					l.attr({'text':dispName});
					break;
				};
				fontSize-=2;
				l.attr({'font-size':fontSize});
			}
		}
		l.data('set',annotationSet)
		//.hover(function(){this.attr({fill:'#F00',cursor:'pointer'});},function(){this.attr({fill:'#000',cursor:'auto'});})
		.hover(function(){setLabelEmphasis(this,'on');},function(){setLabelEmphasis(this,'off');})
		.dblclick(function(){if (confirm('Delete annotation \''+setName+'\'?')) {deleteAnnotation(this)}});

		annotationSet.aspect=l;
		var usedTargetIds={};
		var colorIdx=0; // 0++ -> 1: 0 not used for annotations. Index reset for each annot (PP 29/0516)
		for (var annot in setData) {
			var targetData=setData[annot].split('::'); // targetlist[::color] optional
			var annotColor;
			if (annotNameToColor[type][annot]) { // already seen
				annotColor=annotNameToColor[type][annot];
			}
			else { // first time seen
				if (targetData[1]) { // custom color
					annotColor=targetData[1];
				}
				else { // Choose color from default list
					//var colorIdx=chooseColorIndex(annotColorIndex[type]);
					colorIdx++;
					if (colorIdx>=colorList.length) colorIdx=1;
					annotColor=colorList[colorIdx];
				}
				annotNameToColor[type][annot]=annotColor;
			}
			var targetLabels=targetData[0].split(',');
			for (var i=0; i<targetLabels.length; i++) {
				//if (usedTargetIds[targetLabels[i]]) continue; // ID/label already used by this annotation
				if (matchPattern) {
					/* loose ID match => loop through all label Ids */
					var matchRegExp=new RegExp(matchPattern.replace('###',targetLabels[i]));
					for (var labelID in labelIdList) {
						//if (usedTargetIds[labelID]) continue; // ID/label already used by this annotation
						if (matchRegExp.exec(labelID)) {
							if (usedTargetIds[labelID]) { // ID/label already used by this annotation
								var annotCell=usedTargetIds[labelID];
								annotCell.name+='\n'+annot;
								annotCell.color='#000000';
								annotCell.aspect.attr({fill:'#000000'});
							}
							else {
								var label=labelIdList[labelID];
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
						var annotCell=usedTargetIds[targetLabels[i]];
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
		var maxRank=0;
		for (var sName in annotationList[type]) {maxRank=Math.max(maxRank,annotationList[type][sName].rank);}
		annotationSet.rank=maxRank+1;
		annotationList[type][setName]=annotationSet;
		updateAnnotationLegends();
		return true;
	}
	function addAnnotationCell(type,label,annotationSet,setName,annot,annotStart,cellSize,annotColor) {
		var annotCell=new AnnotationCell(annotationSet,annot,label,annotColor);
		/* Drawing annotation cell */
		var x,y,w,h;
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
		var aC=HM.paper.rect(x,y,w,h,3).attr({fill:annotColor,stroke:'none'}) //,'fill-opacity':0.75 stroke:colorList[colorIdx],'stroke-width':0,'stroke-opacity':0.2
		.data('annotCell',annotCell)
		.hover(function(e){highlightAnnotCell(this,'on',e)},function(){highlightAnnotCell(this,'off')});
		annotCell.aspect=aC;
		label.annotations[setName]=annotCell;
		annotationSet.cells.push(annotCell);
		return annotCell;
	}
	function deleteAnnotation(l) {
		var annotationSet=l.data('set');
		var setName=annotationSet.name;
		var type=annotationSet.type;
		if (annotationOnDelete) {
			var res=annotationOnDelete(type,setName);
			if (!res) {
				alert('An external error has occured. Cannot delete \''+setName+'\'');
				return;
			}
		}

		var annotRank=annotationSet.rank;
		// Delete annotation cells
		for (var i=0; i<annotationSet.cells.length; i++) {
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
		var cellSize=(type=='row')? Math.min(hCellW,20) : (hasCells)? Math.min(hCellH,20) : 20; // 1/2D for type=column (same in AddAnnotation)
		var dx=0,dy=0;
		if (type=='row') {dx=(rowLabelLocation=='left')? -cellSize-2 : cellSize+2;}
		else {dy=(colLabelLocation=='top')? -cellSize-2 : cellSize+2;}
		for (var sName in annotationList[type]) {
			if (annotationList[type][sName].rank > annotRank) {
				annotationList[type][sName].rank--;
				var l2=annotationList[type][sName].aspect;
				if (type=='row') {
					//l2.rotate(-colLabelAngle,l2.attr('x'),l2.attr('y')); // unrotate
					l2.transform(''); // unrotate
					l2.attr({x:l2.attr('x')+dx});
					//l2.rotate(colLabelAngle,l2.attr('x'),l2.attr('y')); // rerotate
					l2.transform('r'+colLabelAngle+','+l2.attr('x')+','+l2.attr('y')); // rerotate
				}
				else {l2.attr({y:l2.attr('y')+dy});}
				for (var i=0; i<annotationList[type][sName].cells.length; i++) {
					var aC=annotationList[type][sName].cells[i].aspect;
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
		var posX=(branchList.row || rowLabelLocation=='right')? paperWidth0 : treeStart.row + 30; //cellX0+mapW+20;
		var posY=(hasCells)? cellY0 : 30; //(hasCells)? cellY0+mapH+20 : 20;
		var types=(hasCells)? ['column','row'] : ['column'];
		for (var t=0; t<2; t++) {
			var type=types[t];
			var numAnnot=0;
			posY+=10;
			for (var setName in annotationList[type]) {
				numAnnot++;
				if (numAnnot==1 && hasCells) {
					annotationLegends.push(HM.paper.text(posX,posY,HM.entities[type]).attr({'font-size':14,'font-weight':'bold','text-anchor':'start'}));
					posY+=15;
				}
				annotationLegends.push(HM.paper.text(posX+10,posY,setName+':').attr({'font-size':12,'font-weight':'bold','text-anchor':'start'}));
				posY+=15;
				var usedCells={};
				for (var i=0; i<annotationList[type][setName].cells.length; i++) {
					var annotName=annotationList[type][setName].cells[i].name;
					if (annotName.match(/\n/)) {annotName='*Multi-label*';}
					else {
						var annotNameInfo=annotName.split('##');
						annotName=annotNameInfo[0];
					}
					if (usedCells[annotName]) {continue;}
					var annotColor=annotationList[type][setName].cells[i].color;
					annotationLegends.push(HM.paper.rect(posX+15,posY-5,10,10).attr({fill:annotColor,stroke:annotColor}));
					annotationLegends.push(HM.paper.text(posX+30,posY,annotName).attr({'font-size':12,'text-anchor':'start'}));
					usedCells[annotName]=1;
					posY+=15;
				}
			}
		}
//console.log(annotationLegends.getBBox());
		/*Readjust paper size if necessary */
		var legend=annotationLegends.getBBox();
        var prevPaperWidth=paperWidth;
        var prevPaperHeight=paperHeight;
		paperWidth=Math.max(paperWidth0,legend.x2+20);
		paperHeight=Math.max(paperHeight0,legend.y2+20);
		HM.paper.setSize(paperWidth,paperHeight);
		var pathArray=backgroundPath.attr('path');
        pathArray[1][1]=paperWidth;
        pathArray[2][1]=paperWidth; pathArray[2][2]=paperHeight;
        pathArray[3][2]=paperHeight;
        backgroundPath.attr({path:pathArray});
//console.log(pathArray);
//console.log(backgroundPath.attr('path'));
		background.attr({width:paperWidth,height:paperHeight});
        var dWidth=paperWidth-prevPaperWidth;
        var dHeight=paperHeight-prevPaperHeight;
        mainDiv.style.width=parseInt(mainDiv.style.width) + dWidth + 'px'; // extend main outer div width
        mainDiv.style.height=parseInt(mainDiv.style.height) + dHeight + 'px'; // extend main outer div height
	}

	function highlightAnnotCell(aC,action,e) {
		var annotCell=aC.data('annotCell');
		if (action=='on') {
			aC.attr({'stroke-width':2,stroke:cellHighlightColor});
			var x=aC.attr('x')+Math.round(aC.attr('width')/2);
			var y=aC.attr('y');
//console.log(aC.transform());
			var textStrg='';
			if (e.ctrlKey || e.shiftKey) {
				textStrg+=HM.entities[annotCell.label.type];
				textStrg+=': '+annotCell.label.label+'\nSet: '+annotCell.set.name+'\n';
			}
			var multiAnnotInfo=annotCell.name.split('\n'); // in case multi annot on same cell
			for (var i=0; i<multiAnnotInfo.length; i++) {
				var annotNameInfo=multiAnnotInfo[i].split('##'); // e.g. 'annotX##1' -> annotX + 1
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
		for (var i=0; i<branchList[type].length; i++) {
			var branch=branchList[type][i];
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
										treePopup[type].hide().transform('')
										treePopup[type][0].attr({'fill-opacity':0.7,'stroke-opacity':0.7}); // box
										treePopup[type][1].attr({'fill-opacity':1}); // text
									}
									);
		}
	}

	function drawTree(type) { // row or column
		if (!branchList[type]) return; // no tree
		var leaves={leaf1:{},leaf2:{}};
		for (var i=0; i<branchList[type].length; i++) {
//console.log(i);
			var branch=branchList[type][i];
			if (type=='row') { // row tree
				for (var leaf in leaves) {
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
				for (var leaf in leaves) {
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
			var branchFullSize=Math.round((branch.distance/treeLength[type])*treeSpace[type]);
			var pathStrg='M'+(leaves.leaf1.x-0.5)+' '+(leaves.leaf1.y-0.5);
			var bNode;
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
			.hover(function(){this.transform('s1.7').attr({'fill-opacity':1,cursor:'pointer'})},function(){if(this.id) {this.attr({'fill-opacity':0,cursor:'auto'}).transform('')}}) //Tree is redrawn if exist sub node => node is deleted during hover!!! (problem for Chrome only so far)
			.click(function(e){if (e.ctrlKey || e.shiftKey) {rotateBranch(this.data('branch'));} else {highLightBranch(this.data('branch'));}});
			//.dblclick(function(){rotateBranch(this.data('branch'));})
			//.click(function(){highLightBranch(this.data('branch'));});
			branch.aspect.push(HM.paper.path(pathStrg).attr({stroke:colorList[branch.selectColorIdx],'stroke-width':1})); // [0] leaves
			branch.aspect.push(bNode); // [1] node
		}
		/* Sending all center circles to front */
		for (var i=0; i<branchList[type].length; i++) {branchList[type][i].aspect[1].toFront();}
		// Tree popups
		if (labelsOnSelect[type]) {
			var t=HM.paper.text(0.5,0.5,labelsOnSelect[type].text).attr({'font-weight':'bold','font-size':14,fill:'#FFF'});
			var b=t.getBBox();
			var tb=HM.paper.rect(b.x-2,b.y-2,b.width+4,b.height+5,2).attr({stroke:'#000',fill:'#000','stroke-opacity':0.7,'fill-opacity':0.7});
			t.toFront();
			treePopup[type]=HM.paper.set(tb,t).hide()
			.hover(function(){this.attr({cursor:'pointer'})},function(){this.attr({cursor:'auto'})})
			.click(function(){focussedLabelsAction('tree',type)}); // former treeClickAction
		}
	}

	function chooseColorIndex(usedColors) {
		var newColorIdx;
		for (var i=1; i<colorList.length; i++) {
			if (!usedColors[i]) {
				newColorIdx=i;
				usedColors[i]=true;
				break;
			}
		}
		if (newColorIdx==null) { // all colors used -> reset
			newColorIdx=1; // 1: already set to true
			//for (var i=2; i<colorList.length; i++) {delete usedColors[i];}
			usedColors={};
			usedColors[1]=true;
		}
		return newColorIdx;
	}

    function highLightBranch(branch,newColorIdx,usedColorIdx) {
		if (newColorIdx==null) {
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
					var translationStrg=(branch.type=='row')? 't'+branch.connectionPos+','+(branch.centerPos-15) : 't'+branch.centerPos+','+(branch.connectionPos-15);
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
		var branchColor=colorList[newColorIdx];
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
		var labelList={leaf1:[],leaf2:[]},numLabels={};
		//List of labels linked to leaf1
		for (var leaf in labelList) {
			if (branch[leaf].centerPos) {pushLabels(branch[leaf],labelList[leaf]);} // points to another branch
			else {labelList[leaf].push(branch[leaf]);} // points to a label
			numLabels[leaf]=labelList[leaf].length;
		}
//console.log(labelList);
		//Swapping leaves
		var cellSize,axis;
		if (branch.type=='row') {
			cellSize=hCellH;
			axis='y';
		}
		else {
			cellSize=hCellW;
			axis='x';
		}
		var shift={leaf1:numLabels.leaf2*cellSize,leaf2:-numLabels.leaf1*cellSize};
		var shiftRanks={leaf1:numLabels.leaf2,leaf2:-numLabels.leaf1};
		var movedLabels=[];
		for (var leaf in shift) {
			for (var i=0; i<labelList[leaf].length; i++) {
				var label=labelList[leaf][i];
				var l=label.aspect;
				//l.attr(axis,l.attr(axis)+shift[leaf]);
				if (branch.type=='row') {l.attr({y:l.attr('y')+shift[leaf]});}
				else {l.attr({x:l.attr('x')+shift[leaf]});} //rotation moved to rank update
				//label.rank+=(leaf=='leaf1')? numLabels.leaf2 : -numLabels.leaf1;
 				movedLabels.push([label,shiftRanks[leaf]]);
				moveCells(label,axis,shift[leaf]);
			}
		}
		//Swapping labels & cells in rowList/columnList data structure
		var movedList,otherList,movedTypeIdx,otherTypeIdx,focusShift=0,focusShift2=0;
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
		var movedCells={};
		for (var i=0; i<movedLabels.length; i++) {
			var label=movedLabels[i][0];
			var rkShift=movedLabels[i][1];
			var oldIdx=label.rank-1;
			label.rank+=rkShift;
			var newIdx=label.rank-1;
			for (var j=0; j<label.heatCellList.length; j++) { // update hCell (row/column)Idx
				label.heatCellList[j][movedTypeIdx]=newIdx;
			}
			movedList[newIdx]=label;
			movedCells[newIdx]=[];
			for (var j=0; j<otherList.length; j++) {
				movedCells[newIdx][j]=otherList[j].heatCellList[oldIdx];
			}
			//Checking if label (& annotations) is inside focusArea
			var trF='',trF2='';
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
			for (var setName in label.annotations) {
				label.annotations[setName].aspect.toFront().transform(trF2);
			}
		}
		for (var idx in movedCells) {
			for (var i=0; i<movedCells[idx].length; i++) {
				otherList[i].heatCellList[idx]=movedCells[idx][i];
				otherList[i].heatCellList[idx][otherTypeIdx]=i;  // update hCell (otherType)Idx
			}
		}
		movedCells=null;

		//Swapping pointers as well
		var newLeaf2=branch.leaf1;
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
		var htmlString='<TABLE cellspacing=0 cellpadding=0><TR>';
		if (HM.editMenu) {
			htmlString+='<TD nowrap><FIELDSET style="padding:2px"><LEGEND><B>Data normalization:</B></LEGEND>';
			/* Scope */
			if (HM.editMenu.scope) {
				htmlString+='<B>Scope:</B><SELECT onchange="cubiojsCharts['+HM.chartID+'].normProcess.scope=this.value; cubiojsCharts['+HM.chartID+'].computeHeatMap()"><OPTION value="row">'+HM.entities.row+'</OPTION><OPTION value="column"';
				if (HM.normProcess.scope=='column') {htmlString+=' selected';}
				htmlString+='>'+HM.entities.column+'</OPTION><OPTION value="all"';
				if (HM.normProcess.scope=='all') {htmlString+=' selected';}
				htmlString+='>Both</OPTION></SELECT>';
			}
			else if (HM.editMenu.scope==false && HM.normProcess.reference != 'z-score') {htmlString+='<B>Scope: '+HM.normProcess.scope+'</B>';}
			/* Reference */
			if (HM.editMenu.type) {
				htmlString+='&nbsp;&nbsp;<B>Reference:</B><SELECT onchange="cubiojsCharts['+HM.chartID+'].normProcess.reference=this.value; cubiojsCharts['+HM.chartID+'].computeHeatMap()"><OPTION value="z-score">z-score</OPTION><OPTION value="mean"';
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
				htmlString+='&nbsp;&nbsp;<B>Color:</B><SELECT onchange="cubiojsCharts['+HM.chartID+'].normProcess.colors=this.value; cubiojsCharts['+HM.chartID+'].computeCellsColor()"><OPTION value="GBR">green/black/red</OPTION><OPTION value="GYR"';
				if (HM.normProcess.colors=='GYR') {htmlString+=' selected';}
				htmlString+='>green/yellow/red</OPTION><OPTION value="GWR"';
				if (HM.normProcess.colors=='GWR') {htmlString+=' selected';}
				htmlString+='>green/white/red</OPTION><OPTION value="BBR"';
				if (HM.normProcess.colors=='BBR') {htmlString+=' selected';}
				htmlString+='>blue/black/red</OPTION><OPTION value="BYR"';
				if (HM.normProcess.colors=='BYR') {htmlString+=' selected';}
				htmlString+='>blue/yellow/red</OPTION><OPTION value="BWR"';
				if (HM.normProcess.colors=='BBR') {htmlString+=' selected';}
				htmlString+='>blue/white/red</OPTION></SELECT>';
			}
			htmlString+='</FIELDSET></TD>\n';
			/* Cell flag visibility */
			if (existFlaggedCells) {
				htmlString+='<TD nowrap><FIELDSET style="padding:2px"><LEGEND><B>Flag:</B></LEGEND><INPUT type="checkbox" checked onchange="cubiojsCharts['+HM.chartID+'].showCellFlag(this.checked)"/><B>Show '+flagText+'&nbsp;</FIELDSET></TD>\n';
			}
		}

		/* Export image button */
		if (HM.exportAsImage) {
			htmlString+='<TD>';
			if (HM.editMenu) htmlString+="&nbsp;&nbsp;";
			if (HM.exportAsImage[3]) { // image format
				htmlString+='<INPUT type="button" value="'+HM.exportAsImage[0]+'" style="font-weight:bold;font-size:10px" onclick="cubiojsCharts['+HM.chartID+'].mergeAndExportToImage(\''+HM.exportAsImage[3]+'\')"/></TD>\n';
			}
			else {
				htmlString+='<FONT style="font-weight:bold">'+HM.exportAsImage[0]+':</FONT>';
				//PNG
				htmlString+='<INPUT type="button" value="PNG" style="font-weight:bold;font-size:10px" onclick="cubiojsCharts['+HM.chartID+'].mergeAndExportToImage(\'png\')"/>';
				//SVG
				htmlString+='<INPUT type="button" value="SVG" style="font-weight:bold;font-size:10px" onclick="cubiojsCharts['+HM.chartID+'].mergeAndExportToImage(\'svg\')"/></TD>\n';
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
		var tmpMap=this.paper.image(canvas.toDataURL(),cellX0,cellY0,mapW,mapH);
		tmpMap.toBack(); // So cell selection (if any) remains visible
		helpLogo.hide();
		helpPopup.hide();
		exportSVGtoImg(paperDivID,this.exportAsImage[1],this.exportAsImage[2],format); // requires chartLibrary2.js
		tmpMap.remove();
		helpLogo.show();
		if (helpIsVisible) helpPopup.show();
	}

    /***********************************************************/
    /************************ Help text ************************/
    /***********************************************************/
    function showHideHelpText() {
		if (helpIsVisible) {
			helpPopup.animate({'fill-opacity':0,'stroke-opacity':0},300,'linear',function() {helpPopup.hide()});
			//helpPopup.hide();
			helpIsVisible=false;
		}
		else {
			helpPopup.show();
			var anim=Raphael.animation({'fill-opacity':0.7,'stroke-opacity':0.7},300,'linear');
			helpPopup[0].animate(anim); // box
			helpPopup[1].animateWith(helpPopup[0],anim,{'fill-opacity':1},300,'linear');

			helpIsVisible=true;
		}
    }

}

/*
####>Revision history<####
# 2.3.5 Displays both PNG and SVG image export options if format is not specified by user (PP 06/06/18)
# 2.3.4 Updates main DIV size when adding/removing annotations legends (PP 11/11/17)
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
