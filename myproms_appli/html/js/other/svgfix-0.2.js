/*
 * svgfix.js - Javascript SVG parser and renderer on Canvas
 * version: 0.2
 * MIT Licensed
 *
 * Fixes svg graphics:
 * - Adds xmlns:xlink
 * - xlink:href instad of href
 * - trim
 *
 * To be used before calls to canvg.js
 *
 * Ignacio Vazquez (ivazquez@adooxen.com)
 * Requires JQuery <--- no longer (PP 14/02/14)
 * http://phpepe.com/
 */
(function(){
	this.svgfix = function (text) {
		var fixed = text ;
		//fixed = jQuery.trim(fixed);
		var xmlnsOcc=fixed.match(/\sxmlns=/g);
//alert('1 xmlns: '+fixed.match(/\sxmlns=/g).length);
		if (xmlnsOcc && xmlnsOcc.length==2) {fixed = fixed.replace(/\sxmlns=\S+/,'');} // A 2nd occurence is added by IE!!!
//alert('2 xmlns: '+fixed.match(/\sxmlns=/g).length);
		if (fixed.indexOf( 'xmlns:xlink' ) == -1 ) {
			fixed = fixed.replace ('<svg ', '<svg xmlns:xlink="http://www.w3.org/1999/xlink" ');
		}
		//fixed = fixed.replace (' href', ' xlink:href');
		fixed = fixed.replace(/\shref/g,' xlink:href');
		//replace dasharray 0 to dasharray 0,0 tos keep error in numberToArray functin in canvg.js
		fixed = fixed.replace(/stroke-dasharray="0"/g, '');
//console.log(fixed);
		return fixed;
	}
})();

 /*
####>Revision history<####
# 0.2.2 Checks for 2nd occurence of xmlns="http:..." added by IE & no longer calls jQuery (PP 14/02/14)
# 0.2.1 add replace stroke-dasharray by nothing to prevent error when exporting svg (SL 17/01/14)
*/