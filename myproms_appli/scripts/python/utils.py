#!/usr/bin/env python3

################################################################################
# utils.py      1.0.1                                                          #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Define some functions useful in various situations and allow to share them   #
# between modules                                                              #
################################################################################
#
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

"""
A bunch of functions that should be useful in various situations.
Functions are made to be adapted in each particular situation.

Functions defined here include : 
- natural_sort
- remove_redundants
- parse_resultsPep
"""

import re

import numpy as np
import pandas as pd


def natural_sort(l):
    """
    Sort the elements of the list in the natural order (not the ASCII one),
    even if the elements contain upper and lower case letters and figures. 
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
    

def remove_redondants(my_list):
    """
    Return a list with elements of my_list, without the redondants elements.
    It keeps the order of first appearance in the list.
    Elements of my_list must be hashable, to be used in the set.
    If not hashable, find an appropriate equivalent.
    """
    
    my_set = set()
    out_list = []
    
    for elem in my_list:
        if elem not in my_set:
            my_set.add(elem)
            out_list.append(elem)
        else:
            continue

    return out_list


def parse_resultsPep(resultsPep_file, filter_out=False):
    """
    Parse resultsPep.txt to extract XIC values of peptides, 
    Transform values from log2 back to natural intensity and 
    Build the matrix of XICs with proteins/peptides as index and states/replicates as columns.
    """
    keep_cols = ['Condition', 'replicate', 'repTech', 'ProteinID', 'Peptide', 'PeptideId', 'log2Measure']
    if filter_out:
        keep_cols.append('out')
    results  = pd.read_csv(resultsPep_file, sep="\t", header='infer', usecols=keep_cols)

    if filter_out:
        results = results.loc[results['out'].isna()]
        results.drop(columns='out', inplace=True)

    results['XIC'] = 2 ** results['log2Measure']
    results.drop(columns='log2Measure', inplace=True)
    
    xic_matrix = results.set_index(['ProteinID', 'Peptide', 'PeptideId', 'Condition', 'replicate', 'repTech'])
    xic_matrix = xic_matrix['XIC'].unstack(level=['Condition', 'replicate', 'repTech'])
    xic_matrix.sort_index(axis=1, level=['Condition', 'replicate', 'repTech'], inplace=True)
    # Sum to remove PeptideId from index, and get one row per peptide 
    # (only one value per line when PeptideId in index, the rest is NaN, so sum is ok to aggregate)
    xic_matrix = xic_matrix.sum(axis=0, level=['ProteinID', 'Peptide'])
    xic_matrix.replace(0, np.nan, inplace=True)

    return xic_matrix

####>Revision history<####
# 1.0.1 [MODIF] Move parse_resultsPep from computeLFQ.py to utils.py (VL 04/05/20)
# 1.0.0 Created (VL 04/05/2020)
