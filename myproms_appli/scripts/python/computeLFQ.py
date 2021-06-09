#!/usr/bin/env python3

################################################################################
# computeLFQ.py      1.0.8                                                     #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Compute the LFQ intensities of proteins, equivalent to MaxQuant LFQ          #
# See Cox et al., MCP 2014
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
python3 computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 2
- minimum 1 peptide for valid ratios (typically when quantifying PTMs), without filtering
python3 computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 1
- LFQ with stabilization of large ratios with minimum of 2 peptides for valid ratios (or stabilization)
python3 computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 2 -s 
- minimum 2 peptides for valid ratios, with filtering
python3 computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 2 -f
- LFQ with minimum 2 peptides for valid ratios and aggregation of replicates into conditions 
python3 computeLFQ.py -i "/path/to/resultsPep.txt" -o "/path/to/LFQ_out.txt" -m 2 -s 

The output is a long format matrix written in a tsv format, with columns "ProteinID Condition Replicate Value".
The "Replicate" column is present only if values are NOT aggragated (-a flag absent)

Note: Computation complexity is O(N^4) with N the number of samples (because of the least squares computation step).
It may take a long time when used with many samples, but it is possible to parallelize the process by groups of proteins
because the results for a given protein does not depend on the other proteins values.
"""

import sys
import argparse
from itertools import combinations

import numpy as np
import pandas as pd
from itertools import combinations

from utils import parse_resultsPep

def setup():
    """
    Parse arguments from command line.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", required=True, help="""Path to the resultsPep.txt file containing the XIC values (assumed to be in log2)""")
    parser.add_argument("-o", "--out_file", required=True, help="""Path to the desired output file which will contain a matrix of LFQ intensities, with protein IDs (or sites) as rows and states as columns""")
    parser.add_argument("-a", "--aggregate_replicates", action="store_true", help="""Whether to aggregate the LFQ values of replicates per condition at the end of the computation.""")
    parser.add_argument("-m", "--min_pep", type=int, default=2, help="""Minimum number of common peptides between two samples to consider the corresponding protein ratio as valid""")
    parser.add_argument("-s", "--stabilize_large_ratios", action="store_true", help="""Whether to use the feature to stabilize large ratios, in the case of samples that differ on a high number of peptides for some proteins.""")
    parser.add_argument("-f", "--filter_out", action="store_true", help="""To use if you want to filter the peptides flagged as 'out' in the file resultsPep.txt""")
    args = parser.parse_args()

    return args


def aggregate_replicates(xic_matrix):
    """
    Aggregate the technical replicates (taking the mean) if there are both biological and technical replicates.
    Otherwise, keep only the level of the ones with multiple replicates (biological or technical).
    If there are neither biological nor technical replicates, just keep the biological replicates index for consistency.
    Note: The aggregate_replicates flag and function are not linked, this is called whether the flag is set or not.
    """
    xic = None
    bio_rep  = xic_matrix.columns.get_level_values('replicate').unique()
    tech_rep = xic_matrix.columns.get_level_values('repTech').unique()
    
    if len(tech_rep) < 1 or len(bio_rep) < 1:
        raise ValueError("resultsPep.txt file is weirdly constructed, there is either no 'replicate' or no 'repTech' values")
    elif len(tech_rep) == 1:
        xic = xic_matrix.droplevel('repTech', axis='columns')
        xic.rename_axis(columns={'replicate': 'Replicate'}, inplace=True)
    elif len(bio_rep) == 1:
        xic = xic_matrix.droplevel('replicate', axis='columns')
        xic.rename_axis(columns={'repTech': 'Replicate'}, inplace=True)
    elif len(bio_rep) > 1 and len(tech_rep) > 1:
        xic = xic_matrix.T.groupby(level=['Condition', 'replicate']).mean().T
        xic.rename_axis(columns={'replicate': 'Replicate'}, inplace=True)
    else:
        raise ValueError("resultsPep.txt file is weirdly constructed, check if the columns 'replicate' and 'repTech' are set properly")
    return xic


def build_ratios_matrix(xic_aggregated, min_pep=2, stabilize_large_ratios=True):
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
    ratios_matrix = None
    xic = xic_aggregated.copy()

    # Map of sample position to name, to avoid dealing with names and having to parse everything
    # and reset columns names to their index
    sample_map = {idx: value for idx, value in zip(range(xic.shape[1]), xic.columns.values)}
    xic.columns = list(range(xic.shape[1]))

    # Get all the possible sample pairs
    sample_pairs = list(combinations(xic.columns, 2))

    # Building the matrix of peptide XIC ratios
    # (used whether stabilization of large ratios is used or not)
    pep_ratios = pd.DataFrame(index=xic.index)  # (peptides x sample pairs)
    for sample_pair in sample_pairs:
        ratio = f"{sample_pair[1]}/{sample_pair[0]}"
        pep_ratios[ratio] = xic[sample_pair[1]] / xic[sample_pair[0]]

    # Aggregation of the peptide ratios into protein ratios, 
    # taking the median ratio of the peptides corresponding to the protein
    prot_groups = pep_ratios.groupby(level='ProteinID')
    median_ratios = prot_groups.median()  # (proteins x sample pairs)

    # Count the number of peptide ratios used to compute the protein ratio
    common_pep_nb = prot_groups.count()  # (proteins x sample pairs)
    # Proteins with nb of peptides ratios < min_pep
    mask_prot = common_pep_nb < min_pep

    if stabilize_large_ratios:
        # Total peptides number per protein per sample 
        total_pep_nb = xic.count(axis='index', level='ProteinID')  # (proteins x samples)

        # Summed peptides XICs (proteins x samples)
        summed_xic = xic.sum(axis='index', level='ProteinID', skipna=True)
    
        # Build the matrix of summed intensities ratios (proteins x sample pairs)
        summed_xic_ratios = pd.DataFrame(index=summed_xic.index)
        # and compute the inverse proportion of common peptides for each sample pair (proteins x sample pairs)
        inv_common_pep_prop = pd.DataFrame(index=common_pep_nb.index, columns=common_pep_nb.columns)
        
        for sample_pair in sample_pairs:
            ratio = f"{sample_pair[1]}/{sample_pair[0]}"
            summed_xic_ratios[ratio] = summed_xic[sample_pair[1]] / summed_xic[sample_pair[0]]
            inv_common_pep_prop[ratio] = total_pep_nb[[sample_pair[0], sample_pair[1]]].max(axis='columns')
        
        summed_xic_ratios.replace([0.0, np.inf], np.nan, inplace=True)
        
        inv_common_pep_prop = inv_common_pep_prop / common_pep_nb
        inv_common_pep_prop.replace(np.inf, np.nan, inplace=True)

        w = (inv_common_pep_prop - 2.5) / 2.5  # weight of linear interpolation for stabilized ratios

        # Compute the stabilized ratios (see "large ratio stabilization" in the MaxLFQ paper)
        stabilized_ratios = pd.DataFrame(index=median_ratios.index, columns=median_ratios.columns)
        stabilized_ratios = median_ratios.where(inv_common_pep_prop < 2.5)
        # Remove protein ratio if nb of peptides ratios < min_pep but only for ratios not stabilized
        stabilized_ratios = stabilized_ratios.mask(mask_prot)
        stabilized_ratios = stabilized_ratios.fillna(summed_xic_ratios.where(inv_common_pep_prop > 5))
        stabilized_ratios = stabilized_ratios.fillna(np.exp(w*np.log(median_ratios) + (1-w)*np.log(summed_xic_ratios)).where((2.5 <= inv_common_pep_prop) & (inv_common_pep_prop <= 5)))

        ratios_matrix = stabilized_ratios.copy()
    else:
        # Remove protein ratio if nb of peptides ratios < min_pep
        ratios_matrix = median_ratios.mask(mask_prot)

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
        # Check if at least one sample has intensity. 
        # If that's the case, sum of proportions = 10 times the number samples with valid ratio(s),
        # otherwise everything equals 0. 
        # 10 times the nb of samples with valid ratios is arbitrary, it's a trade-off between :
        # a value too small which gives very small intensities that are then considered as missing values and
        # a value too big which leads to infinite values at the exponentiation step.
        matrix.iloc[-1, -1] = 10 * ((matrix.iloc[-1, :-1] != 0).sum()) * ((matrix.iloc[-1, :-1] != 0).any().astype(float))

        return matrix
        

    nb_samples = len(sample_map)
    nb_ratios = ratios.shape[1]
    systems_idx = pd.MultiIndex.from_product([ratios.index.values, np.arange(nb_ratios + 1)], names=['ProteinID', 'ratios'])
    systems_col = pd.MultiIndex.from_tuples([sample_map[i] for i in range(nb_samples)] + [('B_matrix', 'B_matrix')], names=['Condition', 'Replicate'])
    
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
    sum_on_prots = sum_on_samples.sum(axis='index', level='ProteinID', skipna=True)

    return sum_on_prots


def compute_LFQ(systems, summed_intensities, aggregate = False):
    """
    For each protein, first compute the solution of the system with least squares.
    Then scale the normalized intensities to the total sum of intensities.
    Finally group the results per state, taking the mean of the replicates, considering null intensities to be missing values (and not as information that the protein is not in the replicate).
    The :aggregate: argument controls whether the replicates are aggregated by condition at the end of the computation.
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
        except np.linalg.LinAlgError:
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

    lfq_intensities = lfq_intensities_sample.copy()
    stack_levels = None
    if aggregate:
        lfq_intensities = lfq_intensities.groupby(axis='columns', level='Condition').mean()
        stack_levels = ["Condition"]
    else:
        stack_levels = ["Condition", "Replicate"]

    lfq_intensities = lfq_intensities.stack(level=stack_levels)
    lfq_intensities.name = "Value"
    lfq_intensities = lfq_intensities.round(decimals=0)

    return lfq_intensities
    

def write_output(lfq_int, output_file):
    lfq_int.to_csv(output_file, sep='\t', na_rep='NA', header=True, index=True)

    return


def main():
    args = setup()

    xic_matrix = parse_resultsPep(args.input_file, args.filter_out)

    xic_aggregated = aggregate_replicates(xic_matrix)

    ratios_matrix, sample_map = build_ratios_matrix(xic_aggregated, args.min_pep, args.stabilize_large_ratios)    

    systems, keep_peptides_from_sample = build_systems_matrix(ratios_matrix, sample_map)

    summed_XIC = sum_peptide_XIC(xic_aggregated, keep_peptides_from_sample)

    lfq_intensities = compute_LFQ(systems, summed_XIC, args.aggregate_replicates)

    write_output(lfq_intensities, args.out_file)

    return


if __name__ == '__main__':
    main()


# Might be good to set 'sort' to False in groupby calls for better performance

####>Revision history<####
# 1.0.8 [FEATURE] Add parameter to choose whether to aggregate replicates into conditions or not (VL 13/01/21)
# 1.0.7 [MINOR] Fix error exception and minor modifs for reading (VL 18/06/20)
# 1.0.6 [BUGFIX] Reimplement usage of min_pep even when using large ratios stabilization (VL 11/06/20)
# 1.0.5 [MODIF] Transfer function for parsing resultsPep.txt to utils.py module to avoid redundancy with computeiBAQ.py (VL 04/06/20)
# 1.0.4 [ENHANCEMENT] Add stabilization of large ratios, adapted from Cox et al. MCP 2014 (VL 29/05/20)
# 1.0.3 [BUGFIX] Adapt the arbitrary value for the least square system to the number of samples (VL 02/03/2020)
# 1.0.2 [MODIF] Compute least-squares solution with (natural) log values instead of non-logged (VL 31/01/20)
# 1.0.1 [BUGFIX] Round LFQ values to the nearest integer, to set very low values of LFQ to 0 (replaced later by NaN), avoid negative values (like -2e-8) and improve readibility (VL 20/01/20)
# 1.0.0 Created (VL 14/01/2020)
