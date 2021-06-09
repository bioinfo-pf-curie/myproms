/*
################################################################################
# idvis.vennDiagram.js    2.0.1                                                #
# uses idvis name space BUT does not require idvis.js itself                   #
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
                  Checks if name space is defined
******************************************************************/
if (typeof(idvis)==='undefined') {
    idvis={};
}


/*****************************************************************
                  Venn diagram object
******************************************************************/
idvis.vennDiagram = function(vennData) {
    //var VD=this;
    const vennDivID=vennData.div;
    var clickAction=vennData.clickAction;
    var proportional=(vennData.proportional===false)? false : true;
    const colors={A:'#DD0000',B:'#0000FF',C:'#C0C000',D:'#00FF00',E:'#555'};
    var usedColors={};
    var vennType=0;
    for (let gr in vennData.groupLabels) {
        vennType++;
        usedColors[gr]=(vennData.groupColors && vennData.groupColors[gr])? vennData.groupColors[gr] : colors[gr];
    }
    var vennW=(!vennData.size)? 400 : (vennData.size < 200)? 200 : vennData.size;
    var vennH=(vennType==2)? Math.round(vennW*0.7) : (vennType==4)? Math.round(vennW*0.8) : vennW;

    /*** Constants ***/
    var fillAttributes={opacity:0.25}; //stroke:'none',
    var strokeAttributes={'stroke-width':2,'stroke-opacity':0.5};
    var labelSize=(vennW >= 450)? 32 : (vennW >= 350)? 24 : 16;
//var labelSpace=labelSize+5;
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
    if (vennData.setValues) {setContents=vennData.setValues;}
    else { // setLists
        setContents={};
    }

    /********************* drawing (Raphael) *********************/
    /* Paper */
    var paper=Raphael(vennDivID,paperW,paperH);
    /* Back panel */
    var bgPanel=paper.rect(0,0,paperW,paperH,0).attr({fill:'#F3F3F3',stroke:'#000','stroke-width':2});
    /* Total */
    var total=0;
    for (let set in setContents) {total+=setContents[set];}
    var grTotal={};
    for (let i=1; i<=vennType; i++) {
        let gr=String.fromCharCode(64+i);
        grTotal[gr]=0;
        for (let set in setContents) {
            if (set.match(gr)) grTotal[gr]+=setContents[set];
        }
    }
    var totTxt=paper.text(2,12,'Total='+total).attr(numAttributes).attr({'text-anchor':'start'});
    if (total===0) {clickAction=null;} // no elements -> no action
    if (clickAction) {
        totTxt.data({label:'all'})
        .hover(function(){highlight(this,'on');},function(){highlight(this,'off');}) // lower-case 'off'
        .click(function() {updateSelection(this);});
    }
    /* Size of a label on paper */
    var lTxt=paper.text(-10,-10,'A').attr(labelAttributes);
    //var labelSpaceX=Math.round(lTxt.getBBox().width/2);
    var labelSpaceY=Math.round(lTxt.getBBox().height/2); // 1/2 because text anchoring in center of x and y
    lTxt.remove();
    var diagRatio=Math.cos(Math.PI/4); // cos(45°)

    /**** 2 sets ****/
    if (vennType==2) {
        if (!setContents.AB) setContents.AB=0;
        let maxR=Math.round(vennW*0.28),
            spacing=0.05*vennW,
            cyA=Math.round(vennH*0.5),cyB=cyA,
            cxA,cxB,rA,rB,dAB,lxA,lyA,lxB,lyB,
            grRatio=grTotal.A/grTotal.B;
        if (grRatio > 20 || grRatio < 1/20) proportional=false;
        if (proportional) {
            /* Circles radius */
            let scale=maxR/Math.sqrt(Math.max(grTotal.A,grTotal.B)/Math.PI); // radius in pixel / radius in gr unit
            rA=Math.sqrt(grTotal.A/Math.PI)*scale;
            rB=Math.sqrt(grTotal.B/Math.PI)*scale;
            /* Circles distance */
            if (setContents.AB) {
                let areaAB=setContents.AB*scale*scale;
                if (grTotal.A > grTotal.B) {
                    if (setContents.B) {dAB=computeGroupDistance(rA,rB,areaAB);}
                    else {dAB=rA-rB;} // no unique to B
                }
                else { // B is bigger
                    if (setContents.A) {dAB=computeGroupDistance(rB,rA,areaAB);}
                    else {dAB=rB-rA;} // no unique to A
                }
            }
            else {dAB=rA+rB+spacing;}
            /* Coordinates on paper & scaling */
            let setWidth=rA+dAB+rB;
            if (setWidth > vennW-2*spacing) {
                let reScale=(vennW-2*spacing)/setWidth;
                rA*=reScale;
                rB*=reScale;
                dAB*=reScale;
                setWidth*=reScale;
            }
            cxA=Math.round((vennW-setWidth)/2+rA);
            cxB=cxA+dAB;
            if (!setContents.A || setContents.A >= setContents.AB) {lxA=cxA; lyA=cyA-rA+labelSpaceY;} else {lxA=Math.round(cxA-(rA-labelSpaceY)*diagRatio); lyA=Math.round(cyA-(rA-labelSpaceY)*diagRatio);}
            if (!setContents.B || setContents.B >= setContents.AB) {lxB=cxB; lyB=cyB-rB+labelSpaceY;} else {lxB=Math.round(cxB+(rB-labelSpaceY)*diagRatio); lyB=Math.round(cyB-(rB-labelSpaceY)*diagRatio);}
        }
        else { // not proportional
            cxA=Math.round(vennW*0.37);
            cxB=vennW-cxA;
            rA=maxR;
            rB=maxR;
            dAB=cxB-cxA;
            lxA=cxA; lyA=cyA-rA+labelSpaceY;
            lxB=cxB; lyB=cyB-rB+labelSpaceY;
        }
        /* Circles */
        drawCircle('A',cxA,cyA,rA);
        drawCircle('B',cxB,cyB,rB);
        /* Labels */
        drawGroupLabel('A',lxA,lyA);
        drawGroupLabel('B',lxB,lyB);
        /* Values */
        let over=(setContents.AB)? rA+rB-dAB : 0; // overlap size
        if (setContents.A) drawValue(true,'A',Math.round(cxA-rA+(2*rA-over)/2),cyA);
        if (setContents.B) drawValue(true,'B',Math.round(cxB+rB-(2*rB-over)/2),cyB);
        if (setContents.AB) {
            let x=Math.round(cxA+rA-over/2);
            if (over > setContents.AB.toString().length*numSize/2) {drawValue(true,'AB',x,cyA);}
            else {drawValue(false,'AB',x,cyA,x,cyA-(cy*0.7));}
        }
        else if (!proportional) {drawValue(true,'AB',Math.round(cxA+rA-over/2),cyA);}
    }
    /**** 3 sets ****/
    else if (vennType==3) {
        let grComb=['A','B','C','AB','AC','BC','ABC'];
        for (let i=0; i<grComb.length; i++) {
            if (!setContents[grComb[i]]) setContents[grComb[i]]=0;
        }
        /* Coordinates */
        let maxR=Math.round(vennW*0.28),
            spacing=0.05*vennW,
            alpha, beta, // Angles at A and B centers
            cxA,cyA,cxB,cyB,cxC,cyC,rA,rB,rC,dAB,dAC,dBC;

/*-----------TEMP-----------*/
if ((setContents.ABC || (setContents.A && setContents.B && setContents.C)) && (setContents.A + setContents.B + setContents.C) / total < 0.4) {proportional=false;}
//console.log((setContents.A + setContents.B + setContents.C) / total);
//proportional=false;
/*--------------------------*/
        if (proportional) {
            /* Coherence check on ABC (if full set overlap) */
            //let sizeUsedABC=((setContents.AB && setContents.AC && setContents.BC) || setContents.ABC)? Math.max(setContents.ABC,Math.round((setContents.AB+setContents.AC+setContents.BC)/3)) : setContents.ABC,
            //let sizeUsedABC=((setContents.AB && setContents.AC && setContents.BC) || setContents.ABC)? Math.max(setContents.ABC,Math.min(setContents.AB,Math.min(setContents.AC,setContents.BC))) : setContents.ABC,
            let duoSet=[setContents.AB,setContents.AC,setContents.BC].sort(),
            //   sizeUsedABC=((setContents.AB && setContents.AC && setContents.BC) || setContents.ABC)? Math.max(setContents.ABC,Math.round((duoSet[0]+duoSet[1])/2)) : setContents.ABC,
               sizeUsedABC=((setContents.AB && setContents.AC && setContents.BC) || setContents.ABC)? Math.max(setContents.ABC,duoSet[0]) : setContents.ABC,
                grTotalUsed={};
sizeUsedABC=setContents.ABC;
console.log('usedABC',setContents.ABC,sizeUsedABC);
            for (let gr in grTotal) {
                grTotalUsed[gr]=grTotal[gr]-setContents.ABC+sizeUsedABC;
            }
            /* Circles radius */
            let scale=maxR/Math.sqrt(Math.max(grTotalUsed.A,Math.max(grTotalUsed.B,grTotalUsed.C))/Math.PI); // radius in pixel / radius in gr uni
            rA=Math.sqrt(grTotalUsed.A/Math.PI)*scale;
            rB=Math.sqrt(grTotalUsed.B/Math.PI)*scale;
            rC=Math.sqrt(grTotalUsed.C/Math.PI)*scale;
            /* Circles distance */
            // A vs B
            if (setContents.AB || setContents.ABC) { // A-B overlap
                let areaAB=(setContents.AB+sizeUsedABC)*scale*scale;
                if (grTotal.A > grTotal.B) {
                    if (setContents.B) {dAB=computeGroupDistance(rA,rB,areaAB);}
                    else {dAB=rA-rB;} // no unique to B
                }
                else { // B is bigger
                    if (setContents.A) {dAB=computeGroupDistance(rB,rA,areaAB);}
                    else {dAB=rB-rA;} // no unique to A
                }
            }
            else {dAB=rA+rB+spacing;}
            // A vs C
            if (setContents.AC || setContents.ABC) { // A-C overlap
                let areaAC=(setContents.AC+sizeUsedABC)*scale*scale;
                if (grTotal.A > grTotal.C) {
                    if (setContents.C) {dAC=computeGroupDistance(rA,rC,areaAC);}
                    else {dAC=rA-rC;} // no unique to C
                }
                else { // C is bigger
                    if (setContents.A) {dAC=computeGroupDistance(rC,rA,areaAC);}
                    else {dAC=rC-rA;} // no unique to A
                }
            }
            else {dAC=rA+rC+spacing;}
            // B vs C
            if (setContents.BC || setContents.ABC) { // B-C overlap
                let areaBC=(setContents.BC+sizeUsedABC)*scale*scale;
                if (grTotal.B > grTotal.C) {
                    if (setContents.C) {dBC=computeGroupDistance(rB,rC,areaBC);}
                    else {dBC=rB-rC;} // no unique to C
                }
                else { // C is bigger
                    if (setContents.B) {dBC=computeGroupDistance(rC,rB,areaBC);}
                    else {dBC=rC-rB;} // no unique to B
                }
            }
            else {dBC=rB+rC+spacing;}
            /* Coordinates on paper */
            cyA=Math.round(vennH*0.37);
            cyB=cyA;
            cxA=Math.round(vennW*0.37);
            cxB=cxA+dAB;
            if ((setContents.AB && setContents.AC && setContents.BC) || setContents.ABC || setContents.A + setContents.B + setContents.C == total) { // full or null overlap
                /* Angle between A-B and A-C lines */
                alpha=Math.acos(((dAC*dAC)+(dAB*dAB)-(dBC*dBC))/(2*dAC*dAB)) || 0; //10*Math.PI/180; // radians NOT degrees ! deg=rad*180/Math.PI
                beta=Math.acos(((dBC*dBC)+(dAB*dAB)-(dAC*dAC))/(2*dAB*dBC)) || 0;
console.log('alpha',alpha,dAB,dAC,dBC,((dAC*dAC)+(dAB*dAB)-(dBC*dBC))/(2*dAC*dAB));
//console.log('alpha',(alpha*180)/Math.PI);
//console.log(dAB,dBC,dAC,alpha,Math.cos(alpha));
                cxC=cxA + dAC * Math.cos(alpha);
                cyC=cyA + dAC * Math.sin(alpha);
                let setOverlap=true,
                    lxA,lyA,lxB,lyB,lxC,lyC;
                if (setContents.A + setContents.B + setContents.C == total) { // No overlap at all
                    setOverlap=false;
                     lxA=cxA; lyA=cyA-rA+labelSpaceY;
                     lxB=cxB; lyB=cyB-rB+labelSpaceY;
                     lxC=cxC; lyC=cyC-rC+labelSpaceY;
                }
                else { // all sets overlap with eachother
                    if (!setContents.A || setContents.A+setContents.AC >= setContents.AB+setContents.ABC) {lxA=cxA; lyA=cyA-rA+labelSpaceY;} else {lxA=Math.round(cxA-(rA-labelSpaceY)*diagRatio); lyA=Math.round(cyA-(rA-labelSpaceY)*diagRatio);}
                    if (!setContents.B || setContents.B+setContents.BC >= setContents.AB+setContents.ABC) {lxB=cxB; lyB=cyB-rB+labelSpaceY;} else {lxB=Math.round(cxB+(rB-labelSpaceY)*diagRatio); lyB=Math.round(cyB-(rB-labelSpaceY)*diagRatio);}
                    lxC=cxC; lyC=cyC+rC-labelSpaceY;
                }

                /** Circles **/
                drawCircle('A',cxA,cyA,rA);
                drawCircle('B',cxB,cyB,rB);
                drawCircle('C',cxC,cyC,rC);

                /** Labels **/
                drawGroupLabel('A',lxA,lyA);
                drawGroupLabel('B',lxB,lyB);
                drawGroupLabel('C',lxC,lyC);

                /** Values with overlap **/
                if (setOverlap) {
                    /* A,B,C values display */
                    let oAB=rA+rB-dAB,
                        oAC=rA+rC-dAC,
                        oBC=rB+rC-dBC,
                        //vyA=(oAC < rA*0.9)? cyA : Math.round(cyA-rA+1.33*(2*rA-oAC)/2),
                        //vyB=(oBC < rB*0.9)? cyB : Math.round(cyB-rB+1.33*(2*rB-oBC)/2),
                        //vyC=(oAC < rC*0.9 && oBC < rC*0.9)? cyC : Math.round(cyC+rC-(2*rC-oBC)/2);
                        vyA=Math.round(cyA-rA+(2*rA-oAC+labelSpaceY)/2),
                        vyB=Math.round(cyB-rB+(2*rB-oBC+labelSpaceY)/2),
                        vyC=Math.round(cyC+rC-(2*rC-oBC+labelSpaceY)/2);
                    drawValue(true,'A',Math.round(cxA-rA+(2*rA-oAB+labelSpaceY)/2),vyA); // A
                    drawValue(true,'B',Math.round(cxB+rB-(2*rB-oAB-labelSpaceY)/2),vyB); // B
                    drawValue(true,'C',cxC,vyC); // C

                    /* Compute circles intersecting points coordinates (2 points) */
                    let intersect={AB:{},AC:{},BC:{}};
                    // AB
                    let coordAB=computeCircleIntersectionPoints(dAB,rA,rB);
                    intersect.AB.x1=intersect.AB.x2=cxA+coordAB.x;
                    intersect.AB.y1=cyA-coordAB.y;
                    intersect.AB.y2=cyA+coordAB.y;
//paper.circle(intersect.AB.x1,intersect.AB.y1,5).attr({stroke:'#F00'});
//paper.circle(intersect.AB.x2,intersect.AB.y2,10).attr({stroke:'#F00'});
                    // AC
                    let coordAC=computeCircleIntersectionPoints(dAC,rA,rC),
                        rCoordAC1=rotateCoordinates(-alpha,coordAC.x,-coordAC.y),
                        rCoordAC2=rotateCoordinates(-alpha,coordAC.x,coordAC.y);
                    intersect.AC.x1=cxA + rCoordAC1.x;
                    intersect.AC.y1=cyA + rCoordAC1.y;
                    intersect.AC.x2=cxA + rCoordAC2.x;
                    intersect.AC.y2=cyA + rCoordAC2.y;
//paper.circle(intersect.AC.x1,intersect.AC.y1,5).attr({stroke:'#0FF'});
//paper.circle(intersect.AC.x2,intersect.AC.y2,10).attr({stroke:'#0FF'});
                    // BC
                    let coordBC=computeCircleIntersectionPoints(dBC,rB,rC),
                        //rCoordBC1=rotateCoordinates(beta,coordBC.x,-coordBC.y),
                        //rCoordBC2=rotateCoordinates(beta,coordBC.x,coordBC.y);
                        rCoordBC1=rotateCoordinates(Math.PI-beta,coordBC.x,-coordBC.y),
                        rCoordBC2=rotateCoordinates(Math.PI-beta,coordBC.x,coordBC.y);
                    intersect.BC.x1=cxB + rCoordBC1.x;
                    intersect.BC.y1=cyB - rCoordBC1.y;
                    intersect.BC.x2=cxB + rCoordBC2.x;
                    intersect.BC.y2=cyB - rCoordBC2.y;
//paper.circle(intersect.BC.x1,intersect.BC.y1,5).attr({stroke:'#00F'});
//paper.circle(intersect.BC.x2,intersect.BC.y2,10).attr({stroke:'#00F'});
                    /* Overlap values display */
                    drawValueInTriangle('AB',intersect.AB.x1,intersect.AB.y1,intersect.AC.x1,intersect.AC.y1,intersect.BC.x1,intersect.BC.y1,cxC,cyC,rC);
                    drawValueInTriangle('AC',intersect.AC.x2,intersect.AC.y2,intersect.BC.x1,intersect.BC.y1,intersect.AB.x2,intersect.AB.y2,cxB,cyB,rB);
                    drawValueInTriangle('BC',intersect.BC.x2,intersect.BC.y2,intersect.AB.x2,intersect.AB.y2,intersect.AC.x1,intersect.AC.y1,cxA,cyA,rA);
                    drawValueInTriangle('ABC',intersect.BC.x1,intersect.BC.y1,intersect.AB.x2,intersect.AB.y2,intersect.AC.x1,intersect.AC.y1);
                }
                /** Values no overlap **/
                else {
                    drawValue(true,'A',cxA,cyA); // A
                    drawValue(true,'B',cxB,cyB); // B
                    drawValue(true,'C',cxC,cyC); // C
                }
            }
            else if (!setContents.ABC) {
                let isLinear=false,
                    setWidth,
                    reScale=function () { // vennW & spacing are global
                        if (setWidth > vennW-2*spacing) {
                            let scale=(vennW-2*spacing)/setWidth;
                            rA*=scale;
                            rB*=scale;
                            rC*=scale;
                            dAB*=scale;
                            dAC*=scale;
                            dBC*=scale;
                        }
                    };
                if (!setContents.AB) { // no A-B overlap => move C in the middle => A-C-B
                    isLinear=true;
                    cyC=cyA;
                    cxC=cxA;
                    cxA=cxC-dAC;
                    cxB=cxC+dBC;
                    setWidth=rA+dAC+dBC+rB;
                    reScale();
                    cxA=spacing+rA;
                    cxC=cxA+dAC;
                    cxB=cxC+dBC;
//console.log('no AB',cxA,cxC,cxB);
                    let oAC=Math.max(0,rA+rC-dAC), // overlap sizes
                        oBC=Math.max(0,rB+rC-dBC),
                    lxA,lyA,lxB,lyB,lxC,lyC;
                    if (!setContents.A || setContents.A >= setContents.AC) {lxA=cxA; lyA=cyA-rA+labelSpaceY;} else {lxA=Math.round(cxA-rA*0.53); lyA=Math.round(cyA-rA*0.63);}
                    lxC=cxC; lyC=cyC-rC+labelSpaceY;
                    if (!setContents.B || setContents.B >= setContents.BC) {lxB=cxB; lyB=cyB-rB+labelSpaceY;} else {lxB=Math.round(cxB+rB*0.53); lyB=Math.round(cyB-rB*0.63);}
                    /* Circles */
                    drawCircle('A',cxA,cyA,rA);
                    drawCircle('B',cxB,cyB,rB);
                    drawCircle('C',cxC,cyC,rC);
                    /* Labels */
                    drawGroupLabel('A',lxA,lyA);
                    drawGroupLabel('B',lxB,lyB);
                    drawGroupLabel('C',lxC,lyC);
                    /* Values */
                    if (setContents.A) {drawValue(true,'A',Math.round(cxA-rA+(2*rA-oAC)/2),cyA);} // A
                    if (setContents.B) {drawValue(true,'B',Math.round(cxB+rB-(2*rB-oBC)/2),cyB);} // B
                    if (setContents.C) {drawValue(true,'C',Math.round(cxC-rC+oAC+(2*rC-oAC-oBC)/2),cyC);} // C
                    if (setContents.AC) {drawValue(true,'AC',Math.round(cxA+rA-oAC/2),cyA);} // AC
                    if (setContents.BC) {drawValue(true,'BC',Math.round(cxB-rB+oBC/2),cyB);} // BC

                }
                else if (!setContents.AC) { // no A-C overlap => move C right of B =>  A-B-C
                    isLinear=true;
                    cxC=cxB+dBC;
                    cyC=cyB;
                    dAC=cxC-cxA; // recompute dAC
                    setWidth=rA+dAB+dBC+rC;
                    reScale();
                    cxA=spacing+rA;
                    cxB=cxA+dAB;
                    cxC=cxB+dBC;
//console.log('no AC');
                    let oAB=rA+rB-dAB, // overlap sizes [cannot be <0 because case "!setContents.AB" already evaluated]
                        oBC=Math.max(0,rB+rC-dBC),
                        lxA,lyA,lxB,lyB,lxC,lyC;
                    if (!setContents.A || setContents.A >= setContents.AB) {lxA=cxA; lyA=cyA-rA+labelSpaceY;} else {lxA=Math.round(cxA-rA*0.53); lyA=Math.round(cyA-rA*0.63);}
                    lxB=cxB; lyB=cyB-rB+labelSpaceY;
                    if (!setContents.C || setContents.C >= setContents.BC) {lxC=cxC; lyC=cyC-rC+labelSpaceY;} else {lxC=Math.round(cxC+rC*0.53); lyC=Math.round(cyC-rC*0.63);}
                    /* Circles */
                    drawCircle('A',cxA,cyA,rA);
                    drawCircle('B',cxB,cyB,rB);
                    drawCircle('C',cxC,cyC,rC);
                    /* Labels */
                    drawGroupLabel('A',lxA,lyA);
                    drawGroupLabel('B',lxB,lyB);
                    drawGroupLabel('C',lxC,lyC);
                    /* Values */
                    drawValue(true,'A',Math.round(cxA-rA+(2*rA-oAB)/2),cyA); // A
                    drawValue(true,'B',Math.round(cxB-rB+oAB+(2*rB-oAB-oBC)/2),cyB); // B
                    drawValue(true,'C',Math.round(cxC+rC-(2*rC-oBC)/2),cyC); // C
                    if (setContents.AB) {drawValue(true,'AB',Math.round(cxA+rA-oAB/2),cyA);} // AB
                    if (setContents.BC) {drawValue(true,'BC',Math.round(cxC-rC+oBC/2),cyB);} // BC
                }
                else if (!setContents.BC) { // no B-C overlap => move C left of A => C-A-B
                    isLinear=true;
                    cxC=cxA-dAC;
                    cyC=cyA;
                    dBC=cxA-cxC; // recompute dBC
                    setWidth=rC+dAC+dAB+rB;
                    reScale();
                    cxC=spacing+rC;
                    cxA=cxC+dAC;
                    cxB=cxA+dAB;
//console.log('no BC');
                    let oAB=rA+rB-dAB, // [cannot be <0 because case "!setContents.AB" already evaluated]
                        oAC=rA+rC-dAC, // [cannot be <0 because case "!setContents.AB" already evaluated]
                        lxA,lyA,lxB,lyB,lxC,lyC;
                    if (!setContents.C || setContents.C >= setContents.AC) {lxC=cxC; lyC=cyC-rC+labelSpaceY;} else {lxC=Math.round(cxC-rC*0.53); lyC=Math.round(cyC-rC*0.63);}
                    lxA=cxA; lyA=cyA-rA+labelSpaceY;
                    if (!setContents.B || setContents.B >= setContents.AB) {lxB=cxB; lyB=cyB-rB+labelSpaceY;} else {lxB=Math.round(cxB+rB*0.53); lyA=Math.round(cyB-rB*0.63);}
                    /* Circles */
                    drawCircle('A',cxA,cyA,rA);
                    drawCircle('B',cxB,cyB,rB);
                    drawCircle('C',cxC,cyC,rC);
                    /* Labels */
                    drawGroupLabel('A',lxA,lyA);
                    drawGroupLabel('B',lxB,lyB);
                    drawGroupLabel('C',lxC,lyC);
                    /* Values */
                    if (setContents.A) {drawValue(true,'A',Math.round(cxA-rA+oAC + (2*rA-oAC-oAB)/2),cyA);} // A
                    if (setContents.B) {drawValue(true,'B',Math.round(cxB-rB+oAB + (2*rB-oAB)/2),cyB);} // B
                    if (setContents.C) {drawValue(true,'C',Math.round(cxC-rC + (2*rC-oAC)/2),cyC);} // C
                    if (setContents.AB) {drawValue(true,'AB',Math.round(cxA+rA-oAB/2),cyA);} // AB
                    if (setContents.AC) {drawValue(true,'AC',Math.round(cxC+rC-oAC/2),cyA);} // AC
                }
                /* Linear Venn */
                if (isLinear) {
                    // Correct Venn height
                    let dH=vennH - (cyA+maxR+10);
                    paperH-=dH;
                    vennH-=dH;
                    paper.setSize(paperW,paperH);
                    bgPanel.attr({height:paperH});
                }
            }
        }
        else { // not proportional
            cxA=Math.round(vennW*0.37); cyA=Math.round(vennH*0.37); rA=maxR;
            cxB=vennW-cxA; cyB=cyA; rB=maxR;
            cxC=Math.round(vennW*0.5); cyC=vennH-cyA; rC=maxR;

            /* Circles */
            drawCircle('A',cxA,cyA,rA);
            drawCircle('B',cxB,cyB,rB);
            drawCircle('C',cxC,cyC,rC);
            /* Labels */
            drawGroupLabel('A',Math.round(cxA-(rA-labelSpaceY)*diagRatio),Math.round(cyA-(rA-labelSpaceY)*diagRatio));
            drawGroupLabel('B',Math.round(cxB+(rB-labelSpaceY)*diagRatio),Math.round(cyB-(rB-labelSpaceY)*diagRatio));
            drawGroupLabel('C',cxC,Math.round(cyC+rC-labelSpaceY));
            /* Values */
            drawValue(true,'A',Math.round(cxA-rA/2),cyA); // A
            drawValue(true,'B',Math.round(cxB+rB/2),cyB); // B
            drawValue(true,'C',cxC,Math.round(cyC+rC*0.4)); // C
            drawValue(true,'AB',Math.round((cxA+cxB)/2),Math.round(cyA-rA*0.33)); // AB
            drawValue(true,'AC',Math.round(cxA-rA*0.1),Math.round(cyA+rA*0.66)); // AC
            drawValue(true,'BC',Math.round(cxB+rB*0.1),Math.round(cyB+rB*0.66)); // BC
            drawValue(true,'ABC',Math.round((cxA+cxB)/2),Math.round(cyB+rA*0.33)); // ABC
        }
    }
    /**** 4 sets ****/
    else if (vennType==4) {
        /* Coordinates */
        let cxA=Math.round(vennW*0.34), cxB=Math.round(vennW*0.48), cxC=vennW-cxB, cxD=vennW-cxA,
            cyAD=Math.round(vennH*0.6), cyBC=Math.round(vennH*0.5),
            rx=Math.round(vennW*0.35), ry=Math.round(vennW*0.17);
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
        let rx=Math.round(vennW*0.35), ry=Math.round(vennW*0.14),
            cxA=Math.round(vennW*0.47), cyA=Math.round(vennH*0.39),
            cxB=Math.round(vennW*0.63), cyB=Math.round(vennH*0.46),
            cxC=Math.round(vennW*0.61), cyC=Math.round(vennH*0.63),
            cxD=Math.round(vennW*0.44), cyD=Math.round(vennH*0.66),
            cxE=Math.round(vennW*0.35), cyE=Math.round(vennH*0.52);
        /* Ellipse */
        drawEllipse('A',cxA,cyA,rx,ry,90,cxA,cyA-rx*0.85);
        drawEllipse('B',cxB,cyB,rx,ry,-18,cxB+rx*0.8,cyB-rx*0.27);
        drawEllipse('C',cxC,cyC,rx,ry,53,cxC+rx*0.5,cyC+rx*0.7);
        drawEllipse('D',cxD,cyD,rx,ry,-53,cxD-rx*0.5,cyD+rx*0.7);
        drawEllipse('E',cxE,cyE,rx,ry,18,cxE-rx*0.8,cyE-rx*0.27);
        /* Values */
        // 1 group
        drawValue(true,'A',cxA,cyA-rx*0.45);
        drawValue(true,'B',cxB+rx*0.4,cyB-rx*0.15);
        drawValue(true,'C',cxC+rx*0.25,cyC+rx*0.35);
        drawValue(true,'D',cxD-rx*0.25,cyD+rx*0.35);
        drawValue(true,'E',cxE-rx*0.4,cyE-rx*0.15);
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


    /* Legends */
    var legend=paper.set();
    if (legendPosition=='side') {
        let y=15;
        for (let gr in vennData.groupLabels) {
            legend.push(paper.text(vennW+10,y,gr+': '+vennData.groupLabels[gr]+' ('+grTotal[gr]+')').attr(numAttributes).attr({'text-anchor':'start',fill:usedColors[gr]}));
            y+=numSize;
        }
    }
    else {
        let y=vennH+10;
        for (let gr in vennData.groupLabels) {
            legend.push(paper.text(10,y,gr+': '+vennData.groupLabels[gr]).attr(numAttributes).attr({'text-anchor':'start',fill:usedColors[gr]}));
            y+=numSize;
        }
    }
    // Update paperW to fit legend if necessary
    var legendBox=legend.getBBox(),
        legendEndX=legendBox.x + legendBox.width + 10,
        legendEndY=legendBox.y + legendBox.height + 10;
    if (legendEndX > paperW) {
        paperW=legendEndX;
    }
    if (legendEndY > paperH) {
        paperH=legendEndY;
    }
    paper.setSize(paperW,paperH);
    bgPanel.attr({width:paperW,height:paperH});

    var popup,
        selectedSet;


    //var selectedSet=(vennData.selected)?
    /* Computes distance between 2 circles center to match overlap wanted */
    function computeGroupDistance(R,r,tgtArea) {
        let d;
        if (tgtArea) { // there is an overlap between the 2 circles
            d=R+r;
            let dd=0.025*r, area=0, prevArea=0,
                R2=R*R, r2=r*r;
            while (area <= tgtArea) {
                prevArea=area;
                d-=dd;
                let d2=d*d;
                area=R2*Math.acos( (d2-r2+R2) / (2*d*R) ) + r2*(Math.acos( (d2+r2-R2) / (2*d*r) ) ) - 0.5 * Math.sqrt( Math.abs( (2*d*R) + (d2-r2+R2) ) * Math.abs( (2*d*R) - (d2-r2+R2) ) );
//console.log('area',d,area);
            }
            if (area-tgtArea > tgtArea-prevArea) { // keep closest match
                d+=dd;
            }
//console.log('Dist',R+r,Math.round(d),tgtArea,area)
        }
        else {d=R+r+5;}
        return Math.round(d);
    }
    /* Compute x,(+/-)y coordinates of the 2 points of 2 intersecting circles (c1Y=c2Y=0, c1X=0) */
    /* From http://mathworld.wolfram.com/Circle-CircleIntersection.html */
    function computeCircleIntersectionPoints(d,r1,r2,rounded) { // d: distance between centers, r1: radius of 1st circle, r2: radius of 2nd circle
        let R=r1, r=r2, reverse=false; // R: radius of larger circle, r: radius of smaller circle
        if (r2 > r1) {R=r2; r=r1; reverse=true;} // use 2nd circle as reference
        let tempV=d*d - r*r + R*R,
            coord={};
        coord.x=tempV / (2*d);
        coord.y=Math.sqrt(4*d*d*R*R - tempV*tempV) / (2*d);
        if (reverse) {coord.x=d-coord.x;}
        if (rounded) {
            coord.x=Math.round(coord.x);
            coord.y=Math.round(coord.y);
        }
        return coord;
    }
    /* Coordinates rotation (https://fr.wikipedia.org/wiki/Rotation_plane) */
    function rotateCoordinates(alpha,x,y) {
        return {
            x : x*Math.cos(alpha) + y*Math.sin(alpha),
            y : -x*Math.sin(alpha) + y*Math.cos(alpha)
        };
    }
    function drawCircle(gr,cx,cy,r) { //,lx,ly
        paper.circle(cx,cy,r).attr(strokeAttributes).attr({stroke:usedColors[gr]}); // Border
        paper.circle(cx,cy,r).attr(fillAttributes).attr({fill:usedColors[gr]}); // Fill
        //let c=paper.circle(cx,cy,r).attr(fillAttributes).attr({fill:usedColors[gr]}).attr({stroke:usedColors[gr]}).attr(strokeAttributes);
//paper.circle(cx,cy,2);
    }
    function drawGroupLabel(gr,lx,ly) {
        let groupText=paper.text(lx,ly,gr).attr(labelAttributes).attr({fill:usedColors[gr]}) // Label
        .data({'label':gr+':full',color:usedColors[gr]})
        .hover(function(){highlight(this,'on');},function(){highlight(this,'off');});
        if (clickAction) {groupText.click(function() {updateSelection(this);});}
        //return paper.set(c1,c2,groupText);
    }
    function drawEllipse(gr,cx,cy,rx,ry,deg,lx,ly) {
        paper.ellipse(cx,cy,rx,ry).attr(strokeAttributes).attr({stroke:usedColors[gr]}).rotate(deg,cx,cy); // Border
        paper.ellipse(cx,cy,rx,ry).attr(fillAttributes).attr({fill:usedColors[gr]}).rotate(deg,cx,cy); // Fill
        let groupText=paper.text(lx,ly,gr).attr(labelAttributes).attr({fill:usedColors[gr]}) // Label
        .data({'label':gr+':full',color:usedColors[gr]})
        .hover(function(){highlight(this,'on');},function(){highlight(this,'off');}); // lower-case 'off'
        if (clickAction) {groupText.click(function() {updateSelection(this);});}
//paper.circle(cx,cy,3);
    }
    /* Finds the "center" of a triangle */
    function drawValueInTriangle(gr,x1,y1,x2,y2,x3,y3,cx,cy,r) { // convex summit, concav summits 1 & 2 (except ABC: all convex)
        let mid23x=x2 + (x3-x2)/2,
            mid23y=y2 + (y3-y2)/2,
            ratio=(gr=='ABC')? 0.67 : 0.4, // true 2/3 rule for ABC -> center of the triangle
            vx=x1 + ratio*(mid23x-x1),
            vy=y1 + ratio*(mid23y-y1),
            d=Math.sqrt((vx-cx)*(vx-cx) + (vy-cy)*(vy-cy)),
            dMax=Math.sqrt((x1-cx)*(x1-cx) + (y1-cy)*(y1-cy)),
            dOpt=r+(dMax-r)/2;
        if (gr != 'ABC' && d < dOpt) { // value is inside opposite circle -> ABC => move outside
 //paper.circle(vx,vy,5);

 //               let dMax=Math.sqrt((x1-cx)*(x1-cx) + (y1-cy)*(y1-cy)),
 //                   dOpt=r+(dMax-r)/2;
 //               ratio*=(d/dOpt);
 //console.log(gr,d,dMax,dOpt,ratio);
 //               vx=x1 + ratio*(mid23x-x1);
 //               vy=y1 + ratio*(mid23y-y1);

            let sx=0.025*(mid23x-x1),
                sy=0.025*(mid23y-y1);
let count=0;
            while (d <= dOpt) {
                ratio-=0.025; //*=0.95;
                vx-=sx;
                vy-=sy;
                d=Math.sqrt((vx-cx)*(vx-cx) + (vy-cy)*(vy-cy));
if (++count==30) break;
            }
//console.log(gr,dMax/r,dOpt,d,count);
        }
//paper.circle(vx,vy,5);
//console.log('triangle',mid23x,mid23y);
        if (gr=='ABC' || dMax-r > 30) {
            drawValue(true,gr,Math.round(vx),Math.round(vy));
        }
        else { // not enough sapce => draw value outside
            let xo,yo;
            if (gr=='AB') {xo=vx; yo=30;}
            else if (gr=='AC') {xo=vx-30; yo=vy+30;}
            else {xo=vx+30; yo=vy+30;}
            drawValue(false,gr,vx,vy,xo,yo);
        }
    }
    function drawValue(inside,gr,inX,inY,outX,outY) {
        let valueText,valueBg,
            value=(setContents[gr])? setContents[gr] : 0; // undefined value
        if (inside) {
            valueText=paper.text(inX,inY,value).attr(numAttributes);
            let b=valueText.getBBox();
            //if (clickAction) valueBg=paper.rect(b.x,b.y+3,b.width,b.height-6,2).attr({fill:'#FFF',stroke:'none',opacity:0}).data({text:valueText}); // better for hover in IE
            valueBg=paper.rect(b.x,b.y+3,b.width,b.height-6,2).attr({fill:'#FFF',stroke:'none',opacity:0}).data({text:valueText}); // better for hover in IE
            valueText.toFront();
        }
        else if (value !== 0) {
            let dot=paper.circle(inX,inY,dotRadius).attr({fill:'black'}),
                line=paper.path('M'+inX+' '+inY+' L'+outX+' '+outY).attr({'stroke-width':1}); // ' L'+x+' '+y
            valueText=paper.text(outX,outY,value).attr(numAttributes);
            let b=valueText.getBBox();
            //valueBg=paper.rect(b.x,b.y+3,b.width,b.height-6,2).attr({fill:'#F3F3F3',stroke:'none'}).data({text:valueText}); // to hide path to dot
            valueBg=paper.rect(b.x,b.y+2,b.width,b.height-4,2).attr({fill:'#F3F3F3',stroke:'none'}).data({text:valueText}); // to hide path to dot
            valueText.toFront()
            .data({dot:dot,line:line});
        }

        if (inside || value > 0) { //clickAction && ()
            valueText.data({label:gr})
            .hover(function(){highlight(this,'on');},function(){highlight(this,'off');}); // lower-case 'off'
            valueBg.hover(function(){highlight(this.data('text'),'on');},function(){highlight(this.data('text'),'off');});
            //if (value > 0) {valueText.click(function() {updateSelection(this)});}
            if (clickAction && value > 0) {valueText.click(function() {updateSelection(this);});}
        }
    }
    function highlight(set,status) {
        let label=set.data('label'),
            color,cursorShape;
        if (status=='on') {
            color='#F00';
            cursorShape='pointer';
            if (label != 'all') {
                let text;
                if (label.match(':')) {
                    let labelInfo=label.split(':');
                    text='All in '+labelInfo[0]+' ('+grTotal[labelInfo[0]]+')';
                }
                else if (label.length==1) {text='Unique to '+label;}
                else {
                    let labelInfo=label.split('');
                    text='Common to '+labelInfo.join(',');
                }
                let popupText=paper.text(set.attr('x'),set.attr('y')-numSize*1.3,text).attr({'font-size':12,'font-weight':'bold',fill:'#F00'}),
                    b=popupText.getBBox(),
                    frame=paper.rect(b.x-2,b.y-1,b.width+4,b.height+2,2).attr({stroke:'#F00',fill:'#FFF','fill-opacity':0.5}),
                    shiftX=(frame.attr('x') < 0)? -1*frame.attr('x')+2 : (frame.attr('x')+frame.attr('width') > paperW)? paperW-(frame.attr('x')+frame.attr('width'))-2 : 0;
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
};

/*
Circle area = Pi*r2, r=sqr(A/Pi)
Area of overlap between 2 circles (From http://calculis.net/comment/aire-intersection-disques-35):
A=R²cos−1( (d² − r² + R²) / 2dR ) + r²( cos−1( (d² + r² − R²) / 2dr ) − ½√( [ 2dR + (d² − r² + R²) ][ 2dR − (d² − r² + R²) ] )
R: radius of biggest circle
r: radius of smaller circle
d: distance between the 2 centres
----
3 sets:
angle in any triangle A-B-C: a²=b²+c²-2bc.cos(alpha) on so on... (a,b,c=opposite segments of A,B,C. alpha=angle at A)
    BC²=AB²+AC²-2*AC*AB*cos(alpha)
    2*AC*AB*cos(alpha)=AB²+AC²-BC²
    alpha=arcos((AC²+AB²-BC²)/(2*AC*AB))
*/

/*
####>Revision history<####
# 2.0.1 [CHANGE] Minor change in color management (PP 07/05/21)
# 2.0.0 Updated vennDiagram.js to code-encapsulated idvis.vennDiagram.js and proportional for 2-3 sets (PP 18/07/17)
# 1.0.6 Proportional cercles for 2 groups (PP ../11/14)
# 1.0.5 Uses paper as variable name instead of canvas (PP 08/10/14)
# 1.0.4 Auto expends paper to fit legends (PP 14/08/14)
# 1.0.3 Minor bug fix when clickAction is not defined (PP 24/04/14)
# 1.0.2 GPL license (PP 23/09/13)
# 1.0.1 Minor changes for cross-browser display optimization (PP 24/07/12)
# 1.0.0 First stable version (PP 19/07/12)
*/