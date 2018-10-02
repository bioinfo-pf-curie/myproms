################################################################################
# promsOboParser.pm    1.0.2                                                   #
# Authors: P. Poullet, G. Arras, F. Yvon (Institut Curie)                      #
# Contact: myproms@curie.fr                                                    #
# Class that inherits from GO::OntologyProvider::OboParser class with some     #
# method redefinitions to fix "has part" relationship problem with OBO files   #
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
package promsOboParser;

use GO::OntologyProvider::OboParser;
use strict;

our @ISA = qw(GO::OntologyProvider::OboParser);

# Fetching key names from parent #
my $PARENT = $GO::OntologyProvider::OboParser::PACKAGE;

my $kNodes          = $PARENT.'::__nodes';
my $kParent         = $PARENT.'::__parent';
my $kAspect         = $PARENT.'::__aspect';

#-------------#
# Constructor #
#-------------#
sub new{
    my ($class,%parameters) = @_;

    if(%parameters){
	return $class->SUPER::new(%parameters);
    } else {
	my $this = {};
	bless($this, $class);
	return $this;
    }
}

#----------------------------------------------------------#
# Redefinitions of GO::OntologyProvider::OboParser methods #
# Theses methods ignore relationships between terms coming #
# from different aspects (P,C,F). This avoids blocking     #
# errors encountered with OboParser in these cases.        #
#----------------------------------------------------------#
sub __populatePaths {

    my $self = shift;

    # go through each GO node in the $kParent hash, the keys of which
    # are the goids that are parents of a given node.

    foreach my $childGoid ( keys %{$self->{$kParent}} ) {

	# note, we directly access the kNodes hash here, rather than
	# use nodeFromId().  This is for performance reasons only -
	# accessing the kNodes hash directly in this method, and the
	# __findAncestor method shaces about 40% of the runtime off of
	# the time taken to populate all the paths.

	my $childNode =  $self->{$kNodes}{$childGoid};

	# now go through each of this child's parents

	#foreach my $parentGoid (@{$self->{$kParent}{$childGoid}}) {

        my $i = 0;
        while($i <= $#{$self->{$kParent}{$childGoid}}){

            my $parentGoid = $self->{$kParent}{$childGoid}[$i];

	    ### Note, there has been a case in the obo file where
	    ### there was an error, and a node was listed as having
	    ### parent in a different aspect.  This results in a fatal
	    ### run time error, as when the parser reads the file, it
	    ### only keeps nodes of a given aspect, and is thus left
	    ### with a dangling reference.  In this case, parentNode
	    ### will be undef, and the call to addParentNodes ends up
	    ### in a run time error.  We can add some logic here to
	    ### give a better error message.

	    my $parentNode = $self->{$kNodes}{$parentGoid}

	    || do {

                # Redefined part #
                # If parent goid has no defined node, then skip it and remove it from parent goid list of the current go term
                splice @{$self->{$kParent}{$childGoid}}, $i, 1;
		warn "There is an error in the obo file, where the relationship between ".
		$childNode->goid.
		" and one or more of its parents is not correctly defined.\n".
		"Please check the obo file.\n".
                next;

	    };

	    ### create connections between child node and its parent

	    $childNode->addParentNodes($parentNode);

	    $parentNode->addChildNodes($childNode);

	    # begin to build the ancestor path, starting with this
	    # parent

	    my @path = ($parentNode);

	    if (exists $self->{$kParent}{$parentGoid}){

		# if this parent has parents, then we continue to
		# build the path upwards to the root.  We pass in the
		# child node, so that each path which reaches the root
		# can be added during the recursive calls to find
		# ancestor

		$self->__findAncestors($childNode,
				       $parentGoid,
				       \@path);

	    }else{

		# otherwise, the path only contains the root, and we add it.

		$childNode->addPathToRoot(@path);

	    }
            $i++;
	}

    }

}

sub __findAncestors {
# Usage:
#
#    $self->__findAncestor($childNode,
#                          $parentGoid,
#                          $pathArrayRef);
#
# This method looks through each goid in hash %{$self->{$kParent}} to
# find all ancestors and push everything to @{$pathArrayRef}..And if
# there is no ancestor found for the $parentGoid, it just add the path
# to the child node.

    my ($self, $childNode, $parentGoid, $pathArrayRef) = @_;

    # go through each immediate parent of the passed in parent

    foreach my $ancestorGoid (@{$self->{$kParent}{$parentGoid}}) {

	# add the ancestor node to our path to the root which is being
	# built

        ### Redefined part ###
        # Skip if ancestor is not a defined node (eg: GO term from an different aspect)
        unless($self->{$kNodes}{$ancestorGoid}){
            next;
        }
        ### -------------- ###

	push (@{$pathArrayRef}, $self->{$kNodes}{$ancestorGoid});

	if (exists $self->{$kParent}{$ancestorGoid}){

	    # if this ancestor has parents, continue building the
	    # paths to the root recursively up the DAG

	    $self->__findAncestors($childNode,
				   $ancestorGoid,
				   $pathArrayRef);

	}else {

	    $childNode->addPathToRoot(reverse @{$pathArrayRef});

	}

	# because there are multiple paths to the root for most nodes,
	# we have now remove the current ancestor from this time
	# through the loop so that the path is reset to the original
	# condition that it was in when passed in to this method

	pop @{$pathArrayRef};

    }

}

#--------------------#
# Additional methods #
#--------------------#
sub sliceAtLevel{
    # Return a sliced version of current ontology at specified level.
    # Actually the sliced ontology contains all nodes of the specified level
    # and their unique parent is rootNode of full ontology.
    my ($this, $level) = @_;

    my @nodes = $this->getNodesAtLevel($level);

    my $slicedOntology = new promsOboParser();
    $slicedOntology->{$kAspect} = $this->__aspect();

    # Filling ontology with OBO formatted node descriptions
    foreach my $node (@nodes) {
	my @entryLine = (
		    'id:'.$node->goid,
		    'name:'.$node->term
		    );
	$slicedOntology->__processNode(\@entryLine);
    }

    $slicedOntology->__populatePaths;

    return $slicedOntology;
}

sub getNodesAtLevel{
    # Return an array of nodes at the specified level of the ontology
    my ($this, $level) = @_;

    my %nodes;

    $this->__getGOatLevel($this->rootNode, $level, \%nodes);

    return values %nodes;
}

sub __getGOatLevel{
    my ($this, $node, $level, $nodeHashRef) = @_;

    if($level > 0 and !$node->isLeaf){
        foreach my $child ($node->childNodes){
            $this->__getGOatLevel($child, $level-1, $nodeHashRef);
        }
    } else {
        $nodeHashRef->{$node->goid} = $node;
    }
}

1;

####>Revision history<####
# 1.0.2 GPL license (PP 23/09/13)
# 1.0.1 New method for slicing ontology at a specified level (FY 02/01/13)
# 1.0.0 New class that inherits from GO::OntologyProvider::OboParser class with some method redefinitions (FY 19/07/12)