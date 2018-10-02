/*
################################################################################
# vennDiagram.js    1.0.5                                                      #
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
                  Venn diagram object
******************************************************************/
function VennDiagram(vennData) {
    var VD=this;
    var vennDivID=vennData.div;
    var clickAction=vennData.clickAction;
    var vennType=0;
    for (var gr in vennData.groupLabels) {
        vennType++;
    }
    var vennW=(!vennData.size)? 400 : (vennData.size < 200)? 200 : vennData.size;
    var vennH=(vennType==2)? Math.round(vennW*0.7) : (vennType==4)? Math.round(vennW*0.8) : vennW;

    /*** Constants ***/
    var colors={A:'#DD0000',B:'#0000FF',C:'#C0C000',D:'#00FF00',E:'#555'};
    var fillAttributes={opacity:0.25}; //stroke:'none',
    var strokeAttributes={'stroke-width':2,'stroke-opacity':0.5};
    var labelSize=(vennW >= 450)? 32 : (vennW >= 350)? 24 : 16;
    var labelAttributes={'font-size':labelSize,'font-weight':'bold'};
    var numSize=(vennW >= 450)? 20 : (vennW >= 350)? 16 : 12;
    var numAttributes={'font-size':numSize,'font-weight':'bold'};
    var dotRadius=(vennW >= 350)? 3 : 2;

    /*** Paper true size ***/
    var paperW,paperH;
    var legendPosition=(vennData.legendPosition && vennData.legendPosition.match('^(side|bottom)$'))? vennData.legendPosition : 'side';
    if (legendPosition=='side') {
        paperW=vennW + 10; // temporary: will be reset to fit legends
        paperH=vennH;
    }
    else {
        paperW=vennW;
        paperH=vennH + 10 + numSize * vennType;
    }

    /*** Sets contents ***/
    var setContents;
    if (vennData.setValues) {setContents=vennData.setValues}
    else { // setLists
        setContents={};
    }

    /********************* drawing (Raphael) *********************/
    /* Paper */
    var paper=Raphael(vennDivID,paperW,paperH);
    /* Back panel */
    var bgPanel=paper.rect(0,0,paperW,paperH,0).attr({fill:'#F3F3F3',stroke:'#000','stroke-width':2});
    /* Legends */
    var legend=paper.set();
    if (legendPosition=='side') {
        var y=15;
        for (var gr in vennData.groupLabels) {
            legend.push(paper.text(vennW+10,y,gr+': '+vennData.groupLabels[gr]).attr(numAttributes).attr({'text-anchor':'start',fill:colors[gr]}));
            y+=numSize;
        }
    }
    else {
        var y=vennH+10;
        for (var gr in vennData.groupLabels) {
            legend.push(paper.text(10,y,gr+': '+vennData.groupLabels[gr]).attr(numAttributes).attr({'text-anchor':'start',fill:colors[gr]}));
            y+=numSize;
        }
    }
    // Update paperW to fit legend if necessary
    var legendBox=legend.getBBox();
    var legendEnd=legendBox.x + legendBox.width + 10;
    if (legendEnd > paperW) {
        paperW=legendEnd;
        paper.setSize(paperW,paperH);
        bgPanel.attr({width:paperW});
    }

    /* Total */
    var total=0;
    for (var set in setContents) {total+=setContents[set];}
    var totTxt=paper.text(2,12,'Total='+total).attr(numAttributes).attr({'text-anchor':'start'});
    if (total==0) {clickAction=null;} // no elements -> no action
    if (clickAction) {
        totTxt.data({label:'all'})
        .hover(function(){highlight(this,'on');},function(){highlight(this,'off')}) // lower-case 'off'
        .click(function() {updateSelection(this)});
    }
    var popup;
    /* Sets */
    var selectedSet;
    /**** 2 sets ****/
    if (vennType==2) {
        /* Coordinates */
        var cxA=Math.round(vennW*0.37), cxB=vennW-cxA;
        var cy=Math.round(vennH*0.5);
        var r=Math.round(vennW*0.28);
        /* Circles */
        drawCircle('A',cxA,cy,r,cxA,Math.round(cy-r*0.8)); // gr,cx,cy,r,lx,ly
        drawCircle('B',cxB,cy,r,cxB,Math.round(cy-r*0.8));
        /* Values */
        drawValue(true,'A',Math.round(cxA-r/2),cy);
        drawValue(true,'B',Math.round(cxB+r/2),cy);
        drawValue(true,'AB',Math.round((cxA+cxB)/2),cy);
    }
    /**** 3 sets ****/
    else if (vennType==3) {
        /* Coordinates */
        var cxA=Math.round(vennW*0.37), cxB=vennW-cxA, cxC=Math.round(vennW*0.5);
        var cyAB=Math.round(vennH*0.37), cyC=vennH-cyAB;
        var r=Math.round(vennW*0.28);
        /* Circles */
        drawCircle('A',cxA,cyAB,r,Math.round(cxA-r*0.53),Math.round(cyAB-r*0.63)); // gr,cx,cy,r,lx,ly
        drawCircle('B',cxB,cyAB,r,Math.round(cxB+r*0.53),Math.round(cyAB-r*0.63));
        drawCircle('C',cxC,cyC,r,cxC,Math.round(cyC+r*0.8));
        /* Values */
        drawValue(true,'A',Math.round(cxA-r/2),cyAB); // A
        drawValue(true,'B',Math.round(cxB+r/2),cyAB); // B
        drawValue(true,'C',cxC,Math.round(cyC+r*0.4)); // C
        drawValue(true,'AB',Math.round((cxA+cxB)/2),Math.round(cyAB-r*0.33)); // AB
        drawValue(true,'AC',Math.round(cxA-r*0.1),Math.round(cyAB+r*0.66)); // AC
        drawValue(true,'BC',Math.round(cxB+r*0.1),Math.round(cyAB+r*0.66)); // BC
        drawValue(true,'ABC',Math.round((cxA+cxB)/2),Math.round(cyAB+r*0.33)); // ABC
    }
    /**** 4 sets ****/
    else if (vennType==4) {
        /* Coordinates */
        var cxA=Math.round(vennW*0.34), cxB=Math.round(vennW*0.48), cxC=vennW-cxB, cxD=vennW-cxA;
        var cyAD=Math.round(vennH*0.6), cyBC=Math.round(vennH*0.5);
        var rx=Math.round(vennW*0.35), ry=Math.round(vennW*0.17);
        /* Ellipse */
        drawEllipse('A',cxA,cyAD,rx,ry,45,Math.round(cxA-rx*0.6),Math.round(cyAD-rx*0.6));
        drawEllipse('B',cxB,cyBC,rx,ry,45,Math.round(cxB-rx*0.55),Math.round(cyBC-rx*0.65));
        drawEllipse('C',cxC,cyBC,rx,ry,-45,Math.round(cxC+rx*0.55),Math.round(cyBC-rx*0.65));
        drawEllipse('D',cxD,cyAD,rx,ry,-45,Math.round(cxD+rx*0.6),Math.round(cyAD-rx*0.6));
        /* Values */
        drawValue(true,'A',Math.round(cxA-rx*0.45),Math.round(cyAD-rx*0.2)); // A
        drawValue(true,'B',Math.round(cxB-rx*0.3),Math.round(cyBC-rx*0.52)); // B
        drawValue(true,'C',Math.round(cxC+rx*0.3),Math.round(cyBC-rx*0.52)); // C
        drawValue(true,'D',Math.round(cxD+rx*0.45),Math.round(cyAD-rx*0.2)); // D
        drawValue(true,'AB',Math.round(cxA-rx*0.18),Math.round(cyAD-rx*0.5)); // AB
        drawValue(true,'AC',Math.round(cxA-rx*0.05),Math.round(cyAD+rx*0.3)); // AC
        drawValue(true,'AD',Math.round((cxA+cxD)/2),Math.round(cyAD+rx*0.62)); // AD
        drawValue(true,'BC',Math.round((cxB+cxC)/2),Math.round(cyBC-rx*0.3)); // BC
        drawValue(true,'BD',Math.round(cxD+rx*0.05),Math.round(cyAD+rx*0.3)); // BD
        drawValue(true,'CD',Math.round(cxD+rx*0.18),Math.round(cyAD-rx*0.5)); // CD
        drawValue(true,'ABC',Math.round(cxA+rx*0.15),Math.round(cyAD-rx*0.2)); // ABC
        drawValue(true,'ABD',Math.round(cxD-rx*0.25),Math.round(cyAD+rx*0.41)); // ABD
        drawValue(true,'ACD',Math.round(cxA+rx*0.25),Math.round(cyAD+rx*0.41)); // ACD
        drawValue(true,'BCD',Math.round(cxD-rx*0.15),Math.round(cyAD-rx*0.2)); // BCD
        drawValue(true,'ABCD',Math.round((cxA+cxD)/2),Math.round(cyAD+rx*0.15)); // ABCD
    }
    /**** 5 sets ****/
    else if (vennType==5) {
//paper.image('/test/Venn5.png',0,0,400,400);
        /* Coordinates */
        var rx=Math.round(vennW*0.35), ry=Math.round(vennW*0.14);
        var cxA=Math.round(vennW*0.47), cyA=Math.round(vennH*0.39);
        var cxB=Math.round(vennW*0.63), cyB=Math.round(vennH*0.46);
        var cxC=Math.round(vennW*0.61), cyC=Math.round(vennH*0.63);
        var cxD=Math.round(vennW*0.44), cyD=Math.round(vennH*0.66);
        var cxE=Math.round(vennW*0.35), cyE=Math.round(vennH*0.52);
        /* Ellipse */
        drawEllipse('A',cxA,cyA,rx,ry,90,cxA,cyA-rx*0.85);
        drawEllipse('B',cxB,cyB,rx,ry,-18,cxB+rx*.8,cyB-rx*0.27);
        drawEllipse('C',cxC,cyC,rx,ry,53,cxC+rx*0.5,cyC+rx*0.7);
        drawEllipse('D',cxD,cyD,rx,ry,-53,cxD-rx*0.5,cyD+rx*0.7);
        drawEllipse('E',cxE,cyE,rx,ry,18,cxE-rx*.8,cyE-rx*0.27);
        /* Values */
        // 1 group
        drawValue(true,'A',cxA,cyA-rx*0.45);
        drawValue(true,'B',cxB+rx*.4,cyB-rx*0.15);
        drawValue(true,'C',cxC+rx*0.25,cyC+rx*0.35);
        drawValue(true,'D',cxD-rx*0.25,cyD+rx*0.35);
        drawValue(true,'E',cxE-rx*.4,cyE-rx*0.15);
        // 2 groups
        drawValue(false,'AB',cxA+rx*0.28,cyA-rx*0.11,cxA+rx*0.5,cyA-rx*0.7); // #1
        drawValue(false,'AC',cxA-rx*0.12,cyA-rx*0.08,cxA*0.65,cyA-rx*0.9); // #5
        drawValue(false,'AD',cxA-rx*0.02,cyA+rx*0.93,cxA-rx*0.02,cyA+rx*1.6); // #3
        drawValue(false,'AE',cxA-rx*0.34,cyA+rx*0.07,cxA-rx*0.6,cyA-rx*0.5); // #5
        drawValue(false,'BC',cxB+rx*0.17,cyB+rx*0.22,cxB+rx*0.83,cyB+rx*0.22); // #2
        drawValue(false,'BD',cxB+rx*0.03,cyB-rx*0.15,cxB+rx*0.25,cyB-rx*0.6); // #1
        drawValue(false,'BE',cxB-rx*0.89,cyB+rx*0.25,cxB-rx*1.52,cyB+rx*0.48); // #4
        drawValue(false,'CD',cxA+rx*0.22,cyA+rx*0.95,cxA+rx*0.22,cyA+rx*1.3); // #3
        drawValue(false,'CE',cxB+rx*0.08,cyB+rx*0.5,cxB+rx*0.8,cyB+rx*0.5); // #2
        drawValue(false,'DE',cxB-rx*0.8,cyB+rx*0.52,cxB-rx*1.55,cyB+rx*0.81); // #4
        // 3 groups
        drawValue(false,'ABC',cxA+rx*0.05,cyA,cxA+rx*0.4,cyA-rx*0.9); // #1
        drawValue(false,'ABD',cxA+rx*0.35,cyA,cxA+rx*0.575,cyA-rx*0.55); // #1
        drawValue(false,'ABE',cxA-rx*0.32,cyA+rx*0.3,cxA-rx*0.6,cyA-rx*0.3); // #5
        drawValue(false,'ACD',cxA+rx*0.12,cyA+rx*0.9,cxA+rx*0.12,cyA+rx*1.4); // #3
        drawValue(false,'ACE',cxA-rx*0.23,cyA+rx*0.04,cxA-rx*0.63,cyA-rx*0.8); // #5
        drawValue(false,'ADE',cxB-rx*0.6,cyB+rx*0.55,cxB-rx*1.6,cyB+rx*0.95); // #4
        drawValue(false,'BCD',cxB+rx*0.03,cyB+rx*0.075,cxB+rx*0.9,cyB+rx*0.075); // #2
        drawValue(false,'BCE',cxB+rx*0.1,cyB+rx*0.35,cxB+rx*0.5,cyB+rx*0.35); // #2
        drawValue(false,'BDE',cxB-rx*0.83,cyB+rx*0.36,cxB-rx*1.4,cyB+rx*0.58); // #4
        drawValue(false,'CDE',cxA+rx*0.325,cyA+rx*0.78,cxA+rx*0.325,cyA+rx*1.6); // #3
        // 4 groups
        drawValue(false,'ABCD',cxA+rx*0.25,cyA+rx*0.1,cxA+rx*0.65,cyA-rx*0.9); // #1
        drawValue(false,'ABCE',cxA-rx*0.12,cyA+rx*0.1,cxA*0.65,cyA-rx*0.7); // #5
        drawValue(false,'ACDE',cxA+rx*0.08,cyA+rx*0.78,cxA+rx*0.08,cyA+rx*1.52); // #3
        drawValue(false,'BCDE',cxB+rx*0,cyB+rx*0.3,cxB+rx*0.7,cyB+rx*0.3); // #2
        drawValue(false,'ABDE',cxB-rx*0.7,cyB+rx*0.4,cxB-rx*1.47,cyB+rx*0.7); // #4
        // 5 groups
        drawValue(true,'ABCDE',cxA+rx*0.05,cyA+rx*0.4);
    }
    // Total

    //var selectedSet=(vennData.selected)?

    function drawCircle(gr,cx,cy,r,lx,ly) {
        paper.circle(cx,cy,r).attr(strokeAttributes).attr({stroke:colors[gr]}); // Border
        paper.circle(cx,cy,r).attr(fillAttributes).attr({fill:colors[gr]}); // Fill
        var groupText=paper.text(lx,ly,gr).attr(labelAttributes).attr({fill:colors[gr]}) // Label
        .data({'label':gr+':full',color:colors[gr]})
        .hover(function(){highlight(this,'on');},function(){highlight(this,'off')});
        if (clickAction) {groupText.click(function() {updateSelection(this)});}
//paper.circle(cx,cy,3);
    }
    function drawEllipse(gr,cx,cy,rx,ry,deg,lx,ly) {
        paper.ellipse(cx,cy,rx,ry).attr(strokeAttributes).attr({stroke:colors[gr]}).rotate(deg,cx,cy); // Border
        paper.ellipse(cx,cy,rx,ry).attr(fillAttributes).attr({fill:colors[gr]}).rotate(deg,cx,cy); // Fill
        var groupText=paper.text(lx,ly,gr).attr(labelAttributes).attr({fill:colors[gr]}) // Label
        .data({'label':gr+':full',color:colors[gr]})
        .hover(function(){highlight(this,'on');},function(){highlight(this,'off')}); // lower-case 'off'
        if (clickAction) {groupText.click(function() {updateSelection(this)});}
//paper.circle(cx,cy,3);
    }
    function drawValue(inside,gr,inX,inY,outX,outY) {
        var valueText,valueBg;
        var value=(setContents[gr])? setContents[gr] : 0; // undefined value
        if (inside) {
            valueText=paper.text(inX,inY,value).attr(numAttributes);
            var b=valueText.getBBox();
            //if (clickAction) valueBg=paper.rect(b.x,b.y+3,b.width,b.height-6,2).attr({fill:'#FFF',stroke:'none',opacity:0}).data({text:valueText}); // better for hover in IE
            valueBg=paper.rect(b.x,b.y+3,b.width,b.height-6,2).attr({fill:'#FFF',stroke:'none',opacity:0}).data({text:valueText}); // better for hover in IE
            valueText.toFront();
        }
        else if (value > 0) {
            var dot=paper.circle(inX,inY,dotRadius).attr({fill:'black'});
            var line=paper.path('M'+inX+' '+inY+' L'+outX+' '+outY).attr({'stroke-width':1}); // ' L'+x+' '+y
            valueText=paper.text(outX,outY,value).attr(numAttributes);
            var b=valueText.getBBox();
            valueBg=paper.rect(b.x,b.y+3,b.width,b.height-6,2).attr({fill:'#F3F3F3',stroke:'none'}).data({text:valueText}); // to hide path to dot
            valueText.toFront()
            .data({dot:dot,line:line})
        }

        if (inside || value > 0) { //clickAction && ()
            valueText.data({label:gr})
            .hover(function(){highlight(this,'on');},function(){highlight(this,'off')}); // lower-case 'off'
            valueBg.hover(function(){highlight(this.data('text'),'on');},function(){highlight(this.data('text'),'off')});
            //if (value > 0) {valueText.click(function() {updateSelection(this)});}
            if (clickAction && value > 0) {valueText.click(function() {updateSelection(this)});}
        }
    }
    function highlight(set,status) {
        var label=set.data('label');
        var color,cursorShape;
        if (status=='on') {
            color='#F00';
            cursorShape='pointer';
            if (label != 'all') {
                var text;
                if (label.match(':')) {
                    var labelInfo=label.split(':');
                    text='All in '+labelInfo[0];
                }
                else if (label.length==1) {text='Unique to '+label;}
                else {
                    var labelInfo=label.split('');
                    text='Common to '+labelInfo.join(',');
                }
                var popupText=paper.text(set.attr('x'),set.attr('y')-numSize*1.3,text).attr({'font-size':12,'font-weight':'bold',fill:'#F00'});
                var b=popupText.getBBox();
                var frame=paper.rect(b.x-2,b.y-1,b.width+4,b.height+2,2).attr({stroke:'#F00',fill:'#FFF','fill-opacity':0.5});
                var shiftX=(frame.attr('x') < 0)? -1*frame.attr('x')+2 : (frame.attr('x')+frame.attr('width') > paperW)? paperW-(frame.attr('x')+frame.attr('width'))-2 : 0;
                if (shiftX) {
                    frame.attr({x:frame.attr('x')+shiftX});
                    popupText.attr({x:popupText.attr('x')+shiftX});
                }
                popupText.toFront();
                popup=paper.set(popupText,frame);
            }
        }
        else { // off or OFF
            if (popup) popup.remove();
            popup=null;
            if (selectedSet && label==selectedSet.data('label') && status=='off') return; //keep highlight if same as selected unless forced (OFF)
            color=(set.data('color'))? set.data('color') : '#000';
            cursorShape='auto';
        }
        set.attr({fill:color,cursor:cursorShape});
        if (set.data('dot')) {
            set.data('dot').attr({fill:color,stroke:color});
            set.data('line').attr({stroke:color});
        }
    }
    function updateSelection(set) {
        if (selectedSet) {
            if (selectedSet.data('label')==set.data('label')) return; // same set
            highlight(selectedSet,'OFF'); // upper-case 'OFF' => force highlight to off even if still hovering
        }
        selectedSet=set;
        clickAction(set.data('label'));
    }
}

/*
####>Revision history<####
# 1.0.5 Uses paper as variable name instead of canvas (PP 08/10/14)
# 1.0.4 Auto expends paper to fit legends (PP 14/08/14)
# 1.0.3 Minor bug fix when clickAction is not defined (PP 24/04/14)
# 1.0.2 GPL license (PP 23/09/13)
# 1.0.1 Minor changes for cross-browser display optimization (PP 24/07/12)
# 1.0.0 First stable version (PP 19/07/12)
*/
