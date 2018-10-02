#!/usr/local/bin/perl -w
$|=1;       # buffer flush (allows instant printing)
################################################################################
# exportSVG.cgi    1.0.2                                                       #
# Authors: P. Poullet, G. Arras, S. Liva (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Export SVG charts to PNG/SVG image for download                              #
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
use strict;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use CGI ':standard';
use MIME::Base64;


my $expName = param('name');
my $imgData = param('imgdata');


if ($expName=~/\.svg/) {
	print header(-type => "image/svg", -attachment => "$expName");
	print $imgData;
}

else { # assume png
	print header(-type => "image/png", -attachment => "$expName");

	#capture, replace any spaces w/ pluses, and decode
	$imgData=~ s/ /+/g;
	my $imgDecoded = MIME::Base64::decode($imgData);
	print $imgDecoded;
}

####> Revision history
# 1.0.2 Export also to svg format (PP 25/06/14)
# 1.0.1 GPL license (PP 11/02/14)