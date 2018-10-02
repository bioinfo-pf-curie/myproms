/*
################################################################################
# barPlot.js    1.0.1                                                          #
# Authors: F. Yvon (Institut Curie)                                            #
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

/*****************************************************
	      barPlot object
 *****************************************************/
function barPlot(param){

      // Graph properties
      this.fixedHeight = (param['height'])?param['height']:300;
      this.fixedWidth = (param['width'])?param['width']:600;
      this.maxWidth = (param['maxWidth'])?param['maxWidth']:0;
      this.minWidth = (param['minWidth'])?param['minWidth']:0;
      this.paddingX = (param['paddingX'])?param['paddingX']:60;
      this.paddingY = (param['paddingY'])?param['paddingY']:35;
      this.barWidth = (param['barWidth'])?param['barWidth']:20;
      this.barColor = (param['barColor'])?param['barColor']:"#0000EE";
      this.onBarClick = (param['onBarClick'])?param['onBarClick']:null;

      // Plotting method
      /* Arguments:
       *	- table that contains key/values following this:
       *		- [key] (id of each data) = {
       *			- ['value'] = quantitative value of data
       *			- ['label'] = text that labels the corresponding bar
       *			- ['color'] = color of the corresponding bar (optional, default is barColor used in constructor)
       *		}
       *	- canvas: the html div id where data should be plotted
       */
      this.draw = function(table, canvas){
	var plot = this;
	var barWidth = plot.barWidth;
	var paddingX = plot.paddingX;
	var paddingY = plot.paddingY;
	var gWidth;
	var gHeight = plot.fixedHeight - 4 * paddingY; // leave 1 padding space at top and 3 at bottom for labels
	var labelPathWidth = 6;
	var onClick = plot.onBarClick;

	// Setting some parameters depending to dataset
	var tLength = size(table);
	var maxValue = max(table);
	var factor = gHeight / maxValue;

	if(this.maxWidth){
	      gWidth = Math.min(tLength*barWidth, plot.maxWidth);
	} else {
	      gWidth = plot.fixedWidth - 2 * paddingX;
	}
	if(this.minWidth){
	      gWidth = Math.max(gWidth, plot.minWidth);
	}

	var paper = Raphael(canvas, gWidth+2*paddingX, plot.fixedHeight);

	// Y and X axis
	var axis = paper.path("M"+paddingX+","+paddingY+"l0,"+gHeight+"l"+gWidth+",0");

	// Y axis labels
	var step;
	if(maxValue <= 10){
	      step = 1;
	} else if (maxValue <= 20){
	      step = 2;
	} else if (maxValue <= 50){
	      step = 5;
	} else {
	      step = Math.pow(10,Math.floor(log10(maxValue)));
	      if(maxValue / step < 3){
		step /= 4;
	      } else if (maxValue / step < 5){
		step /= 2;
	      }
	}
	var originY = gHeight + paddingY;
	for(var i=0;i<=maxValue;i+=step){
	      currentY = originY - i * factor;
	      paper.path("M"+paddingX+","+currentY+"l-"+labelPathWidth+",0");
	      var text = paper.text(paddingX/2, currentY, i);
	}


	// Bar plotting
	var squareWidth = gWidth/tLength;
	if(barWidth > squareWidth){
	  barWidth = squareWidth;
	}
	var currentX = paddingX;
	for(var key in table){
	      var value = table[key]['value'] * factor;
	      var centerX = currentX + (squareWidth/2);
	      var color = (table[key]['color'])?table[key]['color']:plot.barColor;
	      var bar = paper.rect(centerX-(barWidth/2), paddingY+gHeight-value, barWidth, value).attr('fill', color);
	      bar.data('valuePopup',paper.text(centerX, paddingY/2 + gHeight - value, table[key]['value']).hide());

	      // Value displaying at bar top
	      var shortName = table[key]['label'].substring(0,10);
	      if(table[key]['label'].length > 10){
		shortName += '...';
	      }
	      var labelText = paper.text(centerX,
			         originY+(paddingY/2),
			         shortName).attr('text-anchor','end').transform('r-45,'+centerX+','+(originY+(paddingY/2))).toBack();
	      bar.data('label', labelText);

	      // Label displaying at bar bottom
	      var labelPopUpText = paper.text(centerX, paddingY+gHeight+paddingY/2, table[key]['label']).toFront().hide();
	      var popupPadding = 4;
	      var labelPopUp = paper.rect(labelPopUpText.getBBox().x-popupPadding,
			          labelPopUpText.getBBox().y-popupPadding,
			          labelPopUpText.getBBox().width+popupPadding*2,
			          labelPopUpText.getBBox().height+popupPadding*2,
			          2).attr('fill', '#FFFFFF').insertBefore(labelPopUpText).hide();
	      labelPopUp.data('text', labelPopUpText);
	      if(labelPopUp.attr('x')+labelPopUp.attr('width') > gWidth + 2 * paddingX){
		labelPopUp.movePopUp(gWidth+2*paddingX - labelPopUp.attr('width'), labelPopUp.attr('y'));
	      } else if (labelPopUp.attr('x') < 0){
		labelPopUp.movePopUp(0, labelPopUp.attr('y'));
	      }
	      bar.data('labelPopup', labelPopUp);
	      bar.data('id', key);

	      bar.hover(function(){ showPopUp(this);},
		    function(){ hidePopUp(this);});

	      if(onClick){
		bar.click(function(){onClick(this.data('id'))});
	      }
	      currentX+=squareWidth;
	}
      }
}
/******************************************************
	      Static subroutines
 *****************************************************/
function showPopUp(bar){
      bar.data('valuePopup').show();
      bar.data('labelPopup').show();
      bar.data('labelPopup').data('text').show();
}
function hidePopUp(bar){
      bar.data('valuePopup').hide();
      bar.data('labelPopup').hide();
      bar.data('labelPopup').data('text').hide();
}
Raphael.el.movePopUp = function(x, y){
      this.attr('x', x);
      this.attr('y', y);
      this.data('text').attr('x', this.attr('x') + this.attr('width')/2);
      this.data('text').attr('y', this.attr('y') + this.attr('height')/2);
}
function size(table){
      var size = 0;
      for(var key in table){
	size++;
      }
      return size;
}
function max(table){
      var max;
      for(var key in table){
	if(max){
	      if(max*1 < table[key]['value']*1){
		max = table[key]['value'];
	      } else {
	      }
	} else {
	      max = table[key]['value'];
	}
      }
      return max;
}
function log10(n){
      return (Math.log(n) / Math.log(10));
}

/*
####>Revision history<####
# 1.0.1 GPL license (PP 23/09/13)
# 1.0.0 new Barplot class using Raphaël library (FY 31/10/12)
*/
