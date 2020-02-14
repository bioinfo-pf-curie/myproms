#!/usr/bin/env python3

################################################################################
# computeLFQ.py      1.0.2                                                     #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Compute the LFQ intensities of proteins, equivalent to MaxQuant LFQ          #
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
Compute the LFQ intensities of proteins for quantification.
The algorithm is equivalent to MaxLFQ, the LFQ computation from MaxQuant, 
except for the 'Delayed Normalization' of fractions.
As input, it requires at least :
    - the path to the resultsPep.txt file created by AnalysisDiffLimma.R 
    - the path to the desired output file
Optional parameters are :
    - the minimum number of peptides common to both samples needed to consider the ratio as valid
    - whether to filter the peptides marked as 'out' in the resultsPep.txt file

Usage examples : 
- LFQ with minimum of 2 peptides for valid ratios, without filtering
python3 test_computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 2
- minimum 1 peptide for valid ratios (typically when quantifying PTMs), without filtering
python3 test_computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 1
- minimum 2 peptides for valid ratios, with filtering
python3 test_computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 2 -f

The output is a matrix written in the output file, of the shape : 
         State1   State2
prot1   LFQ_1_1  LFQ_1_2
prot2   LFQ_2_1  LFQ_2_2
prot3   LFQ_3_1  LFQ_3_2
prot4   LFQ_4_1  LFQ_4_2

Caution !
In the output matrix, the states are lexicographically ordered, like 
State1, State10, State11, State2, State20, State3
Also, computation complexity is O(N^4) with N the number of samples (because of the least squares computation step)
"""

import sys
import argparse

import numpy as np
import pandas as pd

def setup():
    """
    Parse arguments from command line.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", required=True, help="""Path to the resultsPep.txt file containing the XIC values (assumed to be in log2)""")
    parser.add_argument("-o", "--out_file", required=True, help="""Path to the desired output file which will contain a matrix of LFQ intensities, with protein IDs (or sites) as rows and states as columns""")
    parser.add_argument("-m", "--min_pep", type=int, default=2, help="""Minimum number of common peptides between two samples to consider the corresponding protein ratio as valid""")
    parser.add_argument("-f", "--filter_out", action="store_true", help="""To use if you want to filter the peptides flagged as 'out' in the file resultsPep.txt""")
    args = parser.parse_args()

    return args


def parse_xic(results_pep_file, filter_out=False):
    """
    Parse resultsPep.txt to extract XIC values of peptides, 
    Transform values from log to natural intensity and 
    Build the matrix of XICs with proteins/peptides as index and states/replicates as columns.
    """
    keep_cols = ['Condition', 'replicate', 'repTech', 'ProteinID', 'Peptide', 'PeptideId', 'log2Measure']
    if filter_out:
        keep_cols.append('out')
    results  = pd.read_csv(results_pep_file, sep="\t", header='infer', usecols=keep_cols)

    if filter_out:
        results = results.loc[results['out'].isna()]
        results.drop(columns='out', inplace=True)

    results['XIC'] = 2 ** results['log2Measure']
    results.drop(columns='log2Measure', inplace=True)
    
    xic_matrix = results.set_index(['ProteinID', 'Peptide', 'PeptideId', 'Condition', 'replicate', 'repTech'])
    xic_matrix = xic_matrix['XIC'].unstack(level=['Condition', 'replicate', 'repTech'])
    xic_matrix.sort_index(axis=1, level=['Condition', 'replicate', 'repTech'], inplace=True)
    xic_matrix = xic_matrix.sum(axis=0, level=['ProteinID', 'Peptide'])  # Sum to remove PeptideId from index, and get one row per peptide (only one value per line when PeptideId in index, the rest is NaN, so sum is ok to aggregate)
    xic_matrix.replace(0, np.nan, inplace=True)

    return xic_matrix


def aggregate_replicates(xic_matrix):
    """
    Aggregate the technical replicates (taking the mean) if there are both biological and technical replicates.
    Otherwise, keep only the level of the ones with multiple replicates (biological or technical).
    If there are neither biological nor technical replicates, just keep the biological replicates index for consistency.
    """
    xic = None
    states   = xic_matrix.columns.get_level_values('Condition').unique()
    bio_rep  = xic_matrix.columns.get_level_values('replicate').unique()
    tech_rep = xic_matrix.columns.get_level_values('repTech').unique()
    
    if len(tech_rep) < 1 or len(bio_rep) < 1:
        raise ValueError("resultsPep.txt file is weirdly constructed, there is either no 'replicate' or no 'repTech' values")
    elif len(tech_rep) == 1:
        xic = xic_matrix.droplevel('repTech', axis='columns')
    elif len(bio_rep) == 1:
        xic = xic_matrix.droplevel('replicate', axis='columns')
        xic.rename_axis(columns={'repTech': 'replicate'}, inplace=True)
    elif len(bio_rep) > 1 and len(tech_rep) > 1:
        xic = xic_matrix.T.groupby(level=['Condition', 'replicate']).mean().T
    else:
        raise ValueError("resultsPep.txt file is weirdly constructed, check if the columns 'replicate' and 'repTech' are set properly")   
    return xic


def build_ratios_matrix(xic, min_pep=2):
    """
    Input is a matrix with XICs from peptides, with pairs of 
    protein and peptides as rows and states/replicates as columns
    From that it builds another matrix containing the XIC ratios 
    between all state and bio_rep/tech_rep combinations, for each peptide.
    At this point it is only peptide XIC ratios. 
    Then it aggregates the peptide ratios into protein ratios between samples, 
    by taking the median of the peptide ratios, 
    under the condition that there are at least :min_pep: peptides in common between the samples.
    """
    # Building the matrix of peptide XIC ratios
    pep_ratios = pd.DataFrame(index=xic.index)
    # Map of sample position to name, to avoid dealing with names and having to parse everything
    sample_map = {list(xic.columns.values).index(value): value for value in xic.columns.values}  
    for i in range(len(xic.columns)):
        for j in range(i+1, len(xic.columns)):
            sample_pair = "{0}/{1}".format(j, i)
            pep_ratios[sample_pair] = xic.iloc[:, j] / xic.iloc[:, i]
    
    # Aggregation of the peptide ratios into protein ratios, 
    # taking the median ratio of the peptides corresponding to the protein
    prot_groups = pep_ratios.groupby(level='ProteinID')
    ratios_matrix = prot_groups.median()
    # Count the number of peptide ratios used to compute the protein ratio
    # and remove protein ratio if nb of peptides ratios < min_pep
    pep_nb_matrix = prot_groups.count()
    mask_prot = pep_nb_matrix < min_pep
    ratios_matrix.mask(mask_prot, inplace=True)
    
    return ratios_matrix, sample_map


def build_systems_matrix(ratios, sample_map):
    """
    From matrix of ratios between samples, build the matrices corresponding to the system to solve with the least squares method, for each protein.
    Returns a single pd.DataFrame with system matrices for each protein concatenated.
    """

    def build_system(prot_matrix, ratios_row):
        """
        Build matrix A of system equation for each prot : AX = B, 
        with r(b/a) = Ib / Ia -> ln(r(b/a)) = ln(Ib) - ln(Ia)
        
        For 4 samples a, b, c, d : X = [Ia, Ib, Ic, Id]

              a   b   c   d
        A = [-1   1   0   0]    B = [ln(r(b/a))]    -ln(Ia) + ln(Ib) = ln(r(b/a))
            [-1   0   1   0]        [ln(r(c/a))]    -ln(Ia) + ln(Ic) = ln(r(c/a))
            [-1   0   0   1]        [ln(r(d/a))]    -ln(Ia) + ln(Id) = ln(r(d/a))
            [ 0  -1   1   0]        [ln(r(c/b))]    -ln(Ib) + ln(Ic) = ln(r(c/b))
            [ 0  -1   0   1]        [ln(r(d/b))]    -ln(Ib) + ln(Id) = ln(r(d/b))
            [ 0   0  -1   1]        [ln(r(d/c))]    -ln(Ic) + ln(Id) = ln(r(d/c))
            [ 1   1   1   1]        [       100]    ln(Ia) + ln(Ib) + ln(Ic) + ln(Id) = 100 (arbitrary)

        If we don't have any valid ratio for a given sample, its coefficients are set to 0, which gives a log of 0 and an intensity of 1, which will be replaced by NaN later.
        """
        matrix = prot_matrix.copy()
        for label, ratio in ratios_row.iteritems():
            if pd.notna(ratio):
                idx = ratios_row.index.get_loc(label)
                states = label.split('/')  # ratios columns are like "1/0 2/0 2/1" which corresponds to the sample pairs
                state2 = int(states[0])
                state1 = int(states[1])
                
                matrix.iloc[idx, state1] = -1
                matrix.iloc[idx, state2] = 1
                matrix.iloc[idx, -1] = np.log(ratio)

        # Fill last row of the matrix with 1 if sample has at least one non-null ratio
        matrix.iloc[-1, :-1] = (~((matrix.iloc[:-1, :-1] == 0).all(axis='index'))).astype(float)
        # Check if at least one sample has intensity
        # If that's the case, sum of proportions is 100 (arbitrary), otherwise everything equals 0
        matrix.iloc[-1, -1] = 100 * ((matrix.iloc[-1, :-1] != 0).any().astype(float))

        return matrix
        

    nb_samples = len(sample_map)
    nb_ratios = ratios.shape[1]
    systems_idx = pd.MultiIndex.from_product([ratios.index.values, np.arange(nb_ratios + 1)], names=['ProteinID', 'ratios'])
    systems_col = pd.MultiIndex.from_tuples([sample_map[i] for i in range(nb_samples)] + [('B_matrix', 'B_matrix')], names=['Condition', 'replicate'])
    
    systems = pd.DataFrame(0.0, index=systems_idx, columns=systems_col)
    
    prot_groups = systems.groupby(level='ProteinID')
    systems = prot_groups.apply(lambda x: build_system(x, ratios.loc[x.index.get_level_values('ProteinID')[0]]))

    # For each protein and sample, should the corresponding peptides XICs be considered in the sum of intensities ?
    new_prot_groups = systems.groupby(level='ProteinID')
    keep_peptides_from_sample = new_prot_groups.apply(lambda x: (x.iloc[-1, :-1].astype(bool)))
    
    return systems, keep_peptides_from_sample


def sum_peptide_XIC(xic_matrix, xics_to_keep):
    """
    For every protein, compute the sum of XICs on all the corresponding peptides and all samples. 
    The summed intensities will later be used to get the protein intensity absolute scale. 
    """
    xic_valid = xic_matrix.mask(~xics_to_keep)
    sum_on_samples = xic_valid.sum(axis='columns', skipna=True)
    sum_on_prots = sum_on_samples.sum(level='ProteinID', skipna=True)

    return sum_on_prots


def compute_LFQ(systems, summed_intensities):
    """
    For each protein, first compute the solution of the system with least squares.
    Then scale the normalized intensities to the total sum of intensities.
    Finally group the results per state, taking the mean of the replicates, considering null intensities to be missing values (and not as information that the protein is not in the replicate).
    """

    def solve_system_with_lstsq(sys_matrix):
        """
        Given the right matrices, solve the system of equations for the log values of sample intensities 
        (normalized to 100) under the least squares constraint. 
        Return the intensity values after exponentiation (de-logging) and rescaling to 1.
        """
        prot_ID = sys_matrix.index.get_level_values('ProteinID')[0]
        A = sys_matrix.drop(columns='B_matrix', level='Condition')
        B = sys_matrix[('B_matrix', 'B_matrix')]

        try:
            X = np.linalg.lstsq(A, B, rcond=None)  # Tuple : (np.ndarray: lstsq solution, np.ndarray: residuals, int: rank of A, np.ndarray: singular values of A)
        except np.LinAlgError as e:
            print("Could not find a solution for the system during computation of LFQ intensities for protein {}.".format(prot_ID), file=sys.stdout)
            print("LFQ intensities will be set to 0 for this protein.", file=sys.stdout)
            X = (                                                   # To be consistent with np.linalg.lstsq output 
                np.zeros(A.shape[1], dtype=float),                  # Least squares solution
                np.zeros(1, dtype=float),                           # residuals
                0,                                                  # matrix rank
                np.zeros(min(A.shape[0], A.shape[1]), dtype=float)  # singular values of A
            )        

        solved_sys_exp = pd.Series(data=np.exp(X[0]), index=A.columns, name=prot_ID)
        # Replace value by 0 if intensity is around 1, because it means log(intensity) was 0, so no valid ratio
        solved_sys_exp.loc[solved_sys_exp.round(decimals=0) == 1.0] = 0.0
        solved_sys_norm = solved_sys_exp / solved_sys_exp.sum()

        return solved_sys_norm

    solved_systems = systems.groupby(level='ProteinID').apply(solve_system_with_lstsq)

    lfq_intensities_sample = solved_systems.multiply(summed_intensities, axis='index')
    lfq_intensities_sample = lfq_intensities_sample.round(decimals=0)  # Avoid values very close to 0 (like -2e-8, which are probably the result of the least square method with floats) and improve readibility
    lfq_intensities_sample.mask(lfq_intensities_sample == 0.0, inplace=True)
    
    lfq_intensities = lfq_intensities_sample.groupby(axis='columns', level='Condition').mean()
    lfq_intensities = lfq_intensities.round(decimals=0)

    return lfq_intensities
    

def write_output(lfq_int, output_file):
    lfq_int.to_csv(output_file, sep='\t', na_rep='NA', header=True, index=True)

    return


def main():
    args = setup()

    xic_matrix = parse_xic(args.input_file, args.filter_out)

    xic_aggregated = aggregate_replicates(xic_matrix)

    ratios_matrix, sample_map = build_ratios_matrix(xic_aggregated, args.min_pep)    

    systems, keep_peptides_from_sample = build_systems_matrix(ratios_matrix, sample_map)

    summed_XIC = sum_peptide_XIC(xic_aggregated, keep_peptides_from_sample)

    lfq_intensities = compute_LFQ(systems, summed_XIC)

    write_output(lfq_intensities, args.out_file)

    return


if __name__ == '__main__':
    main()


# Might be good to set 'sort' to False in groupby calls for better performance

####>Revision history<####
# 1.0.2 [MODIF] Compute least-squares solution with (natural) log values instead of non-logged (VL 31/01/20)
# 1.0.1 [BUGFIX] Round LFQ values to the nearest integer, to set very low values of LFQ to 0 (replaced later by NaN), avoid negative values (like -2e-8) and improve readibility (VL 20/01/20)
# 1.0.0 Created (VL 14/01/2020)