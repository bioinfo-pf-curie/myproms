#!/usr/bin/env python3

################################################################################
# proteomic_ruler.py      1.0.0                                                #
# Authors: V. Laigle (Institut Curie)                                          #
# Contact: myproms@curie.fr                                                    #
# Performs absolute quantification of proteins, given intensity values         #
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
Estimate cellular copy numbers and concentrations from protein intensities using the proteomic ruler approach.
See Wisniewski, Hein, et al., MCP, 2014. PMID 25225357 for more information on the Proteomic Ruler
This code is adapted from the proteomic ruler implementation of Perseus (Jurgen Cox lab).
(https://github.com/JurgenCox/perseus-plugins/blob/master/PluginProteomicRuler/CopyNumbers.cs)
"""

import argparse
import math
import os
import re
import sys
from operator import itemgetter

import numpy as np
import pandas as pd

AVOGADRO = 6.02214129e23  # mol^(-1)
BASE_PAIR_WEIGHT = 615.8771  # Da

input_df = None
prot_ids_idx = None
prot_ids = None
intensities_idx = None
intensities = None
logarithmized = None
log_base = None
averaging_mode = None
groups_file = None
molecular_weights = None
detectability_correction = None
correction_factor = None
scaling_mode = None
protein_amount_per_cell = None
ploidy = None
custom_proteins = None
custom_prot_qty = None
protein_concentration = None
organism_name = None
output = None
output_file = None


def setup():
    global input_df
    global prot_ids_idx
    global prot_ids
    global intensities_idx
    global intensities
    global logarithmized
    global log_base
    global averaging_mode
    global groups_file
    global molecular_weights
    global detectability_correction
    global correction_factor
    global scaling_mode
    global protein_amount_per_cell
    global ploidy
    global custom_proteins
    global custom_prot_qty
    global protein_concentration
    global organism_name
    global output
    global output_file

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--num_parameters_file", required=False, default=None, help="""TSV file containing the 
    numerical 
    parameters""")
    parser.add_argument("-c", "--char_parameters_file", required=True, help="""TSV file containing the string 
    parameters""")
    args = parser.parse_args()

    if args.num_parameters_file:
        num_parameters = parse_parameters(args.num_parameters_file)
        char_parameters = parse_parameters(args.char_parameters_file)
        params = Namespace(**num_parameters, **char_parameters)
    else:
        char_parameters = parse_parameters(args.char_parameters_file)
        params = Namespace(**char_parameters)

    input_df = pd.read_csv(params.input_matrix, sep='\t', header='infer')  # pd.DataFrame
    prot_ids_idx = params.protein_acc
    if prot_ids_idx is None or prot_ids_idx < 0:
        print("Invalid protein IDs index, provide the index of the protein IDs column please", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(0)
    prot_ids = input_df.iloc[:, prot_ids_idx]  # pd.Series, because only one column selected

    intensities_idx = params.intensities  # List of int, even if it contains only one index
    if not intensities_idx:
        print("No intensities index, provide some please", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(0)
    intensities = input_df.iloc[:, intensities_idx]  # pd.DataFrame because columns are given in a list

    logarithmized = params.logarithmized
    if logarithmized:
        log_base = params.log_base
        if not log_base:
            print("Missing or incorrect log_base whereas logarithmized columns were indicated.", file=sys.stderr)
            print("Provide the right logarithm used on the intensities (2, natural or 10).", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(0)

    averaging_mode = params.averaging_mode
    if averaging_mode == 2:
        try:
            groups_file = params.groups_file
        except AttributeError:
            print("Missing groups_file whereas the same normalization within groups was asked for.", file=sys.stderr)
            print("Provide the file that gives the samples' groups.", file=sys.stderr)
            print("Exiting...")
            sys.exit(0)
        else:
            if not groups_file or not os.path.isfile(groups_file):
                print("Missing groups_file whereas the same normalization within groups was asked for.", file=sys.stderr)
                print("What was provided is not a file, please provide the file that gives the samples' groups.", file=sys.stderr)
                print("Exiting...")
                sys.exit(0)

    molecular_weights = params.molecular_weights

    detectability_correction = params.detectability_correction
    if detectability_correction:
        try:
            correction_factor = params.correction_factor_idx
        except AttributeError:
            print("Missing correction factor whereas detectability correction was asked for .", file=sys.stderr)
            print("Provide the column index of correction factor.", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(0)
        else:
            if not correction_factor or correction_factor < 0:
                print("Missing or incorrect correction factor whereas detectability correction was asked for .", file=sys.stderr)
                print("Provide the column index of correction factor.", file=sys.stderr)
                print("Exiting...", file=sys.stderr)
                sys.exit(0)

    if params.total_protein_amount:
        scaling_mode = 0
        protein_amount_per_cell = params.total_protein_amount
    elif params.histone_proteomic_ruler:
        scaling_mode = 1
        ploidy = params.histone_proteomic_ruler
    else:
        scaling_mode = 2
        custom_proteins = params.custom_proteins
        custom_prot_qty = params.custom_prot_qty
        if len(custom_prot_qty) != len(intensities_idx):
            print("You didn't provide the custom protein amounts corresponding to the intensity columns, please do so.", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(0)

    protein_concentration = params.protein_concentration
    organism_name = params.organism_name
    output = params.output  # List of int from [0, 1, 2, 3, 4, 5, 6, 7, 8]

    if correction_factor is not None:
        indexes = [prot_ids_idx, molecular_weights, correction_factor] + intensities_idx
    else:
        indexes = [prot_ids_idx, molecular_weights] + intensities_idx
    if any(idx > (len(input_df.columns) - 1) for idx in indexes):
        print("There is a problem with the matrix and indexes provided", file=sys.stderr)
        print("Some indexes are greater than the number of columns of the matrix", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(0)

    output_file = params.out_file

    return


def parse_parameters(parameter_file):
    parameters = {}
    with open(parameter_file, 'r') as param_file:
        for line in param_file:
            key, value = line.strip().split('\t', 1)

            if key == 'input_matrix':  # TSV file with protein IDs, MWs, intensities for each sample etc.
                parameters[key] = value

            elif key == 'protein_acc':  # Protein identifiers (as Uniprot ACC) column index. It is considered as
                # protein IDs in this script but they are not the protein IDs from MyProMS (which may be confusing)
                parameters[key] = int(value)

            elif key == 'intensities':  # semicolon-separated list of intensity columns indexes
                parameters[key] = [int(col) for col in value.strip().split(';')]

            elif key == 'logarithmized':  # Boolean
                if value.lower() in ['1', 'true', 't', 'yes', 'y']:
                    parameters[key] = True
                else:
                    parameters[key] = False

            elif key == 'log_base':  # Required only if logarithmized is True
                if value in ['2', 'natural', '10']:
                    parameters[key] = value
                else:
                    parameters[key] = None

            elif key == 'averaging_mode':
                parameters[key] = int(value)

            elif key == 'groups_file':
                if value.lower() in ['', 'none', 'undef', 'na', 'n/a', 'nan']:
                    parameters[key] = None
                else:
                    parameters[key] = value

            elif key == 'molecular_weights':
                parameters[key] = int(value)

            elif key == 'detectability_correction':  # Boolean
                if value.lower() in ['1', 'true', 't', 'yes', 'y']:
                    parameters[key] = True
                else:
                    parameters[key] = False

            elif key == 'correction_factor_idx':
                if value.lower() in ['', 'none', 'undef', 'na', 'n/a', 'nan']:
                    parameters[key] = None
                else:
                    parameters[key] = int(value)

            elif key == 'total_protein_amount':
                if value.lower() in ['', 'none', 'undef', 'na', 'n/a', 'nan']:
                    parameters[key] = None
                else:
                    parameters[key] = float(value)

            elif key == 'histone_proteomic_ruler':
                if value.lower() in ['', 'none', 'undef', 'na', 'n/a', 'nan']:
                    parameters[key] = None
                else:
                    parameters[key] = float(value)

            elif key == 'custom_proteins':
                if value.lower() in ['', 'none', 'undef', 'na', 'n/a', 'nan']:
                    parameters[key] = None
                else:
                    parameters[key] = value.strip().split(';')

            elif key == 'custom_prot_qty':
                if value.lower() in ['', 'none', 'undef', 'na', 'n/a', 'nan']:
                    parameters[key] = None
                else:
                    parameters[key] = [float(qty) for qty in value.strip().split(';')]

            elif key == 'protein_concentration':
                parameters[key] = float(value)

            elif key == 'organism_name':
                if value.lower() in ['', 'none', 'undef', 'na', 'n/a', 'nan']:
                    parameters[key] = None
                else:
                    parameters[key] = value

            elif key == 'output':
                parameters[key] = [int(out_code) for out_code in value.strip().split(';')]

            elif key == 'out_file':
                parameters[key] = value

            else:
                pass

    return parameters


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return str(self.__dict__)


def process_data():
    output_df = input_df.copy()
    suppl_tables = []

    columns = intensities.copy(deep=True)
    input_names = []
    sample_names = []
    cust_prot_qty = custom_prot_qty

    # Hold all intensity values
    intensity_regex = re.compile(r'^(?:(?:LFQ )?[Ii]ntensity )?(.*?)(?: \[Intensities\])?$')
    for col in intensities:
        input_names.append(col)
        sample_names.append(intensity_regex.search(col).group(1))

    # Average over columns if option selected (and replace previous columns)
    if averaging_mode == 3:
        # Original code (perseus plugin proteomic ruler) considers zeros as valid as in
        # the following (commented) line:
        # columns = intensities.median(axis=1, skipna=True).to_frame()
        # We don't consider zeros as valid here
        # axis=1: columns axis, apply function to each row. Don't consider NA or zeros
        columns = intensities.where(intensities != 0).median(axis=1, skipna=True).to_frame()
        if cust_prot_qty:
            cust_prot_qty = [np.median(cust_prot_qty)]
        sample_names = [""]

    # Revert logarithm if necessary
    if logarithmized:
        def exp2(x):
            return math.pow(2, x)

        def exp10(x):
            return math.pow(10, x)

        log_bases = {'2': exp2, 'natural': math.exp, '10': exp10}
        log_reverse = np.vectorize(log_bases[log_base])
        if 0 in columns.to_numpy():
            print("Are the columns really logarithmized ?", file=sys.stderr)
            print("They contain zeroes!", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(0)
        columns = columns.applymap(log_reverse)

    # Deal with missing values by replacing it with 0 after checking/reversing logarithm
    # No incidence on the following steps
    columns.fillna(0, inplace=True)

    mw = output_df.iloc[:, molecular_weights]  # pd.Series
    # Define whether the molecular masses are given in Da or kDa
    if mw.median(skipna=True) < 250:  # Most likely kDa, convert it to Da
        mw = mw.multiply(1000)

    detectability_norm_factor = mw.copy()  # pd.Series
    if detectability_correction:
        if not isinstance(correction_factor, int):
            print("Missing correction factors column index (or a wrong type of data was provided)", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(0)
        detectability_norm_factor = output_df.iloc[:, correction_factor]  # pd.Series
    # The normalization factor needs to be nonzero for all proteins
    # Check and replace with 1 for all relevant cases
    for row_label, row in detectability_norm_factor.iteritems():
        if row == 0 or pd.isna(row):
            detectability_norm_factor.loc[row_label] = 1

    # Detect the organism with protein IDs or get the name directly from parameters
    if organism_name:
        organism = get_organism_by_name(organism_name)
    else:
        organism = detect_organism(prot_ids.to_list())
    # C-value: amount of DNA per haploid genome, see: http://en.wikipedia.org/wiki/C-value
    c_value = organism.genome_size * BASE_PAIR_WEIGHT / AVOGADRO

    if scaling_mode == 1:
        # Find the histones
        histone_rows = find_histones(prot_ids.to_list(), organism)
        # Write a categorical column indicating the histones
        #histone_cat_col = []

        #for row in range(len(input_df.index)):
        #    histone_cat_col.append("+" if row in histone_rows else "")
        #histone_series = pd.Series(histone_cat_col, index=input_df.index)

        #output_df["Histones"] = histone_series

    # Do the same for custom proteins
    elif scaling_mode == 2:
        custom_protein_rows = find_custom_proteins(prot_ids.to_list(), custom_proteins)
        #custom_prot_cat_col = []

        #for row in range(len(input_df.index)):
        #    custom_prot_cat_col.append("+" if row in custom_protein_rows else "")
        #custom_prot_series = pd.Series(custom_prot_cat_col, index=input_df.index)

        #output_df["Custom proteins"] = custom_prot_series

    normalization_factors = []

    # Calculate normalization factor for each column
    for col_label, col in columns.iteritems():
        col_idx = columns.columns.get_loc(col_label)
        # Normalization factor to go from intensities to copies,
        # needs to be determined using the right scaling approach:
        # total protein, histone proteomic ruler, custom protein approach
        if scaling_mode == 0:  # total protein
            mw_weighted_normalized_intensities = col * mw / detectability_norm_factor
            summed_mw_weighted_normalized_intensities = mw_weighted_normalized_intensities.sum(skipna=True)
            factor = protein_amount_per_cell * 1e-12 * AVOGADRO / summed_mw_weighted_normalized_intensities

        elif scaling_mode == 1:  # histone proteomic ruler
            mw_normalized_histone_intensities = (col.iloc[histone_rows] * mw.iloc[histone_rows] /
                                                 detectability_norm_factor.iloc[histone_rows])
            summed_mw_normalized_histone_intensities = mw_normalized_histone_intensities.sum(skipna=True)
            factor = c_value * ploidy * AVOGADRO / summed_mw_normalized_histone_intensities

        elif scaling_mode == 2:  # custom protein approach
            mw_normalized_custom_protein_intensities = (col.iloc[custom_protein_rows] * mw.iloc[custom_protein_rows] /
                                                        detectability_norm_factor.iloc[custom_protein_rows])
            summed_mw_normalized_custom_protein_intensities = mw_normalized_custom_protein_intensities.sum(skipna=True)
            factor = cust_prot_qty[col_idx] * 1e-12 * AVOGADRO / summed_mw_normalized_custom_protein_intensities

        else:
            print("Couldn't retrieve the scaling mode. Default normalization factor is 1.")
            factor = 1.0

        normalization_factors.append(factor)

    # Check averaging mode
    # averaging_mode == 0 --> All columns separately, nothing to change
    if averaging_mode == 1:  # Same factor for all
        factor = sum(normalization_factors) / float(len(normalization_factors))
        normalization_factors = [factor] * len(normalization_factors)

    if averaging_mode == 2:  # Same factor in each group
        sample_groups = []
        with open(groups_file, 'r') as groups_f:
            for line in groups_f:
                try:
                    sample, group = line.strip().split('\t')
                except ValueError:
                    print("The groups file provided doesn't have the expected format.", file=sys.stderr)
                    print("It should be a TSV file with only two columns: "
                          "samples (column 1) and their respective group (column 2)", file=sys.stderr)
                    print("If a sample is not assigned any group, state NA as the group", file=sys.stderr)
                    print("Exiting...", file=sys.stderr)
                    sys.exit(0)
                else:
                    if sample in sample_names:
                        if (not group) or (group.lower() in ["na", "n/a", "nan"]):
                            sample_groups.append((sample, ""))
                        else:
                            sample_groups.append((sample, group))
                    else:
                        continue

        if len(sample_groups) != len(columns.columns):
            print("The number of samples in the groups' file is different from "
                  "the number of intensity columns provided (or the sample names are different).", file=sys.stderr)
            print("It cannot work, each sample should correspond to exactly one "
                  "intensity column and should be assigned a group.", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(0)

        # Re-order sample_groups by sample names to get the same order as the intensity columns
        names_sample_groups = [elem[0] for elem in sample_groups]
        indexes_sample_groups = [names_sample_groups.index(name) for name in sample_names]
        sample_groups = [sample_groups[idx] for idx in indexes_sample_groups]

        group_names = [elem[1] for elem in sample_groups]
        unique_group_names = unique_in_list(group_names)
        grouping = []
        nb_groups = len(unique_group_names)
        for group in group_names:
            if group and group in unique_group_names:
                grouping.append(unique_group_names.index(group))
            else:  # If invalid group, assign a different value than for the other groups.
                # If the groups are assigned properly, the program should never enter this else statement
                # It will enter if some groups are NA for example
                grouping.append(nb_groups)
                nb_groups += 1

        factors = {}
        for i in range(len(columns.columns)):
            if grouping[i] in factors.keys():
                factors[grouping[i]].append(normalization_factors[i])
            else:
                factors[grouping[i]] = [normalization_factors[i]]

        average_normalization_factors = []
        for i in range(len(columns.columns)):
            factor = factors[grouping[i]]
            average_normalization_factors.append(sum(factor) / float(len(factor)))

        normalization_factors = average_normalization_factors

    # Initialize the variables for the annotation rows
    sample_name_row = pd.Series(None, index=input_df.columns)
    input_name_row = pd.Series(None, index=input_df.columns)
    total_protein_row = pd.Series(None, index=input_df.columns)
    total_molecules_row = pd.Series(None, index=input_df.columns)
    # Populate the organism_row variable with N/A strings as defaults
    # Not None, because pandas considers it as a float,
    # which may cause errors when writing the annotations in the end.
    organism_row = pd.Series("N/A", index=input_df.columns)

    if scaling_mode == 1:
        histone_mass_row = pd.Series(None, index=input_df.columns)
        ploidy_row = pd.Series(None, index=input_df.columns)
    elif scaling_mode == 2:
        custom_prot_mass_row = pd.Series(None, index=input_df.columns)
    cell_volume_row = pd.Series(None, index=input_df.columns)

    # Loop over all selected columns and calculate copy numbers
    for col_label, col in columns.iteritems():
        col_idx = columns.columns.get_loc(col_label)
        sample_name = sample_names[col_idx]
        factor = normalization_factors[col_idx]

        copy_numbers = col / detectability_norm_factor * factor
        total_molecules = copy_numbers.sum(axis=0, skipna=True)
        protein_mass = copy_numbers * mw * 1e12 / AVOGADRO  # picograms
        total_protein = protein_mass.sum(axis=0, skipna=True)  # picograms
        if scaling_mode == 1:
            histone_mass = protein_mass.iloc[histone_rows].sum(axis=0, skipna=True)  # picograms
        elif scaling_mode == 2:
            custom_prot_mass = protein_mass.iloc[custom_protein_rows].sum(axis=0, skipna=True)  # picograms

        total_volume = total_protein / protein_concentration * 1000  # femtoliters

        suffix = "" if sample_name == "" else " {}".format(sample_name)
        if 0 in output:
            column_name = "Copy number{}".format(suffix)
            output_df[column_name] = copy_numbers.round(0)
        if 1 in output:
            concentrations = copy_numbers / (total_volume * 1e-15) / AVOGADRO * 1e9  # nanomolar
            column_name = "Concentration [nM]{}".format(suffix)
            output_df[column_name] = concentrations.round(3)
        if 2 in output:
            mass_per_cell = copy_numbers * mw * 1e12 / AVOGADRO  # picograms
            column_name = "Mass per cell [pg]{}".format(suffix)
            output_df[column_name] = mass_per_cell.round(4)
        if 3 in output:
            mass_fraction = copy_numbers * mw * 1e12 / AVOGADRO / total_protein * 1e6  # ppm
            column_name = "Abundance (mass/total mass) [*10^-6]{}".format(suffix)
            output_df[column_name] = mass_fraction.round(3)
        if 4 in output:
            mole_fraction = copy_numbers / total_molecules * 1e6  # ppm
            column_name = "Abundance (molecules/total molecules) [*10^-6]{}".format(suffix)
            output_df[column_name] = mole_fraction.round(3)

        if 5 in output or 6 in output:
            ranks = copy_numbers.rank(method='min', na_option='keep', ascending=False)
            non_valid_copy_numbers = ((copy_numbers == 0) | (copy_numbers.isna()))
            valid_ranks = ranks[~non_valid_copy_numbers].count()
            ranks[non_valid_copy_numbers] = np.nan

            if 5 in output:
                column_name = "Copy number rank{}".format(suffix)
                output_df[column_name] = ranks.round(0)
            if 6 in output:
                relative_ranks = ranks / float(valid_ranks)
                column_name = "Relative copy number rank{}".format(suffix)
                output_df[column_name] = relative_ranks.round(3)

        if intensities_idx[col_idx] < len(input_df.index) and averaging_mode != 3:
            input_name_row.iloc[intensities_idx[col_idx]] = input_names[col_idx]
            sample_name_row.iloc[intensities_idx[col_idx]] = sample_names[col_idx]
            total_protein_row.iloc[intensities_idx[col_idx]] = round(total_protein, 2)
            total_molecules_row.iloc[intensities_idx[col_idx]] = round(total_molecules, 0)
            organism_row[intensities_idx[col_idx]] = organism.name
            if scaling_mode == 1:
                histone_mass_row.iloc[intensities_idx[col_idx]] = round(histone_mass, 4)
                ploidy_row.iloc[intensities_idx[col_idx]] = round(histone_mass * 1e-12 / c_value, 2)
            elif scaling_mode == 2:
                custom_prot_mass_row.iloc[intensities_idx[col_idx]] = round(custom_prot_mass, 4)
            cell_volume_row.iloc[intensities_idx[col_idx]] = round(total_volume, 2)  # femtoliters

    # Summary annotation row
    if averaging_mode != 3 and 7 in output:
        output_df.loc["Total protein [pg/cell]"] = total_protein_row.round(2)
        output_df.loc["Total molecules per cell"] = total_molecules_row.round(0)
        output_df.loc["Organism"] = organism_row
        if scaling_mode == 1:
            output_df.loc["Histone mass [pg/cell]"] = histone_mass_row.round(4)
            output_df.loc["Ploidy"] = ploidy_row.round(2)
        elif scaling_mode == 2:
            output_df.loc["Custom protein mass [pg/cell]"] = custom_prot_mass_row.round(4)
        output_df.loc["Cell volume [fl]"] = cell_volume_row.round(2)

    # Summary matrix
    if averaging_mode != 3 and 8 in output:
        suppl_table = pd.DataFrame({"Sample": sample_name_row,
                                    "Input Column": input_name_row,
                                    "Organism": organism_row,
                                    "Total protein [pg/cell]": total_protein_row.round(2),
                                    "Total molecules per cell": total_molecules_row.round(0),
                                    "Cell volume [fl]": cell_volume_row.round(2)
                                    })
        if scaling_mode == 1:
            suppl_table["Histone mass [pg/cell]"] = histone_mass_row.round(4)
            suppl_table["Ploidy"] = ploidy_row.round(2)
        elif scaling_mode == 2:
            suppl_table["Custom protein mass [pg/cell]"] = custom_prot_mass_row.round(4)

        suppl_tables.append(suppl_table)

    return output_df, suppl_tables


def unique_in_list(my_list):
    results = []
    unique_el = set()
    for elem in my_list:
        if elem in unique_el:
            continue
        else:
            results.append(elem)
            unique_el.add(elem)
    return results


class Organism:
    """
    An object representing a model organism
    """

    def __init__(self, name=None, other_names=None, genome_size=None, histone_ids=None):
        self.name = name if name else "n.d."
        self.other_names = other_names if other_names else []
        self.genome_size = genome_size if genome_size else None
        self.histone_ids = histone_ids if histone_ids else []
        self.custom_prot_ids = []

    def __eq__(self, organism):
        return isinstance(organism, Organism) and self.name == organism.name

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return "Organism: {0}, {1}, {2}".format(self.name, self.genome_size, self.other_names)

    def __repr__(self):
        return "Organism: {0}, {1}, {2}".format(self.name, self.genome_size, self.other_names)

def supported_organisms():
    """
    The list of the organisms that are supported.
    These organisms and their histones can be auto-detected, provided that uniprot IDs are used.
    Returns a list of Organism objects
    """

    organisms = []
    h_sapiens = Organism(
            name="H. sapiens",
            other_names=["human", "h. sapiens", "homo sapiens", "h sapiens", "homo sapiens sapiens",
                         "h sapiens sapiens", "h. sapiens sapiens"],
            genome_size=3200000000,
            histone_ids=["P07305", "Q8IZA3", "Q92522", "P0C5Y9", "P0C5Z0", "H0YFX9", "Q9BTM1", "A8MQC5", "C9J0D1",
                         "C9J386", "E5RJU1", "Q71UI9", "P16104", "B4DJC3", "D6RCF2", "O75367", "Q5SQT3", "Q9P0M6",
                         "P0C0S5", "P0C1H6", "A9UJN3", "P57053", "Q7Z2G1", "B4DEB1", "P84243", "B2R4P9", "K7EMV3",
                         "K7ES00", "K7EK07", "K7EP01", "Q6NXT2", "Q02539", "P16401", "P16403", "P16402", "Q4VB24",
                         "P10412", "A3R0T8", "A1L407", "P22492", "Q96QV6", "P04908", "Q08AJ9", "Q93077", "P20671",
                         "P0C0S8", "A3KPC7", "Q96KK5", "Q99878", "A4FTV9", "Q92646", "Q96A08", "P33778", "P62807",
                         "P58876", "B2R4S9", "Q93079", "P06899", "O60814", "Q99880", "I6L9F7", "Q99879", "Q99877",
                         "P23527", "P68431", "P62805", "Q99525", "Q0VAS5", "B2R4R0", "Q6FI13", "Q8IUE6", "Q16777",
                         "Q16778", "B4DR52", "Q5QNW6", "Q71DI3", "Q5TEC6", "Q7L7L0", "Q8N257", "Q16695", "Q6TXQ4",
                         "Q14463", "B4E0B3", "B2R5B6", "A2RUA4", "B2R5B3", "Q9HA11", "A8K9J7", "B2R6Y1", "B4E380",
                         "A8K4Y7", "Q6B823", "Q6LBZ2", "A3R0T7"]
    )

    organisms.append(h_sapiens)
    m_musculus = Organism(
            name="M. musculus",
            other_names=["mouse", "m. musculus", "mus musculus", "m musculus"],
            genome_size=2700000000,
            histone_ids=["Q9DAD9", "B2RTM0", "Q8CBB6", "Q921L4", "Q5M8Q2", "Q810S6", "B1AV31", "Q497L1", "A9Z055",
                         "Q8CGP9", "P10922", "Q8CJI4", "E0CZ52", "E0CYL2", "Q8VIK3", "Q80ZM5", "Q9CQ70", "Q8R1M2",
                         "Q3THW5", "Q8R029", "B2RVP5", "P27661", "Q9QZQ8", "Q8CA90", "Q8BP16", "Q9CTR1", "Q8CCK0",
                         "Q9D3V6", "Q9D3U7", "Q3UA95", "Q3TFU6", "G3UWL7", "G3UX40", "P0C0S6", "F8WI35", "E0CZ27",
                         "E0CYN1", "E0CYR7", "P84244", "P02301", "Q9QYL0", "P43275", "P43276", "P15864", "Q5SZA3",
                         "P43277", "Q149Z9", "P43274", "Q07133", "I7HFT9", "Q8CGP4", "P22752", "B2RVF0", "Q61668",
                         "Q8CGP5", "A0AUV1", "Q8CGP6", "A3KPD0", "Q8CGP7", "F8WIX8", "A0JNS9", "P70696", "Q64475",
                         "Q6ZWY9", "P10853", "Q64478", "A0JLV3", "Q8CGP1", "B2RVD5", "P10854", "B2RTK3", "Q8CGP2",
                         "P68433", "P84228", "A1L0U3", "A1L0V4", "P62806", "B2RWH3", "Q6GSS7", "Q64522", "Q64523",
                         "Q149V4", "Q64525", "G3X9D5", "Q64524", "B9EI85", "Q61667", "Q8BFU2", "A2AB79", "Q9D2U9",
                         "Q8CGP0", "Q6B822", "P07978", "Q9D9Z7"]
    )
    organisms.append(m_musculus)
    d_melanogaster = Organism(
            name="D. melanogaster",
            other_names=["drosophila", "fruitfly", "fruit fly", "d. melanogaster", "drosophila melanogaster",
                         "d melanogaster"],
            genome_size=130000000,
            histone_ids=["Q6TXQ1", "P02255", "Q4AB54", "Q4ABE3", "Q4ABD8", "Q4AB94", "P84051", "Q4AB57", "P08985",
                         "P02283", "P02299", "E2QCP0", "P84249", "P84040"]
    )
    organisms.append(d_melanogaster)

    c_elegans = Organism(
            name="C. elegans",
            other_names=["worm", "c. elegans", "caenorhabditis elegans", "c elegans", "roundworm", "nematode"],
            genome_size=100300000,
            histone_ids=["P10771", "P15796", "Q19743", "O17536", "O01833", "Q9U3W3", "Q18336", "P09588", "J7S164",
                         "J7SA65", "Q27485", "Q23429", "Q27511", "P04255", "Q27894", "P08898", "K7ZUH9", "Q10453",
                         "Q9U281", "Q27490", "Q27532", "P62784", "Q27484", "Q27876", "O16277", "Q27489"]
    )
    organisms.append(c_elegans)

    s_cerevisiae = Organism(
            name="S. cerevisiae",
            other_names=["budding yeast", "s. cerevisiae", "s cerevisiae", "saccharomyces cerevisiae"],
            genome_size=12100000,
            histone_ids=["P53551", "P04911", "P04912", "Q12692", "P02293", "P02294", "P61830", "P02309"]
    )
    organisms.append(s_cerevisiae)

    s_pombe = Organism(
            name="S. pombe",
            other_names=["fission yeast", "s. pombe", "s pombe", "schizosaccharomyces pombe"],
            genome_size=14100000,
            histone_ids=["P48003", "P04909", "P04910", "P04913", "P09988", "P10651", "P09322"]
    )
    organisms.append(s_pombe)

    g_gallus = Organism(
            name="G. gallus",
            other_names=["chicken", "g. gallus", "g gallus", "gallus gallus", "red junglefowl"],
            genome_size=1230000000,
            histone_ids=["O93327", "P62801", "P84229", "P0C1H3", "P02259", "P0C1H5", "P0C1H4", "Q9PSW9", "P02272",
                         "P02263", "P84247", "P70081", "P15340", "P35062", "P70082", "P08286", "Q6XXM1", "P08284",
                         "P09987", "P08287", "Q5ZMD6", "P08288", "P08285", "P84553"]
    )
    organisms.append(g_gallus)

    return organisms

# function supported_organism_names from original perseus plugin removed (was not used)


def get_organism_by_name(name):
    """
    Maps the organism name to a conventional name and gets the organism object from there 

    Parameter: name of the organism
    Returns the detected organism
    :type name: str
    """

    if "_" in name:
        org_name = name.replace("_", " ").lower()
    else:
        org_name = name.lower()

    organisms = supported_organisms()

    for organism in organisms:
        if org_name in organism.other_names:
            return organism
        else:
            continue
    print("Organism was not recognized, please try something else", file=sys.stderr)
    print("Exiting...", file=sys.stderr)
    sys.exit(0)


def detect_organism(protein_group_ids):
    """
    Finds the organism given a set of Protein IDs
    The function will look at all the protein IDs and match them to
    the list of known histone IDs for each supported organism
    The organism with most hits will be returned
    
    Parameter: pd.Series of protein group IDs (strings of semicolon-separated IDs)
    Returns the detected organism
    """

    histone_hits = {}
    for organism in supported_organisms():
        histone_hits[organism] = 0

    for protein_group_id in protein_group_ids:
        protein_ids = protein_group_id.strip().split(';')
        for prot_id in protein_ids:
            for organism in supported_organisms():
                if prot_id in organism.histone_ids:
                    histone_hits[organism] += 1

    if max(histone_hits.values()) == 0:
        return Organism()
    organisms_counts = sorted(histone_hits.items(), key=itemgetter(1), reverse=True)

    return organisms_counts[0][0]


def find_histones(protein_group_ids, organism):
    """
    Finds those protein groups that represent histones in the given organism
    Parameters:
        - list of protein group IDs (strings of semicolon-separated IDs)
        - the organism (as an Organism object)
    Returns the list of row indices that represent histones
    """

    histone_rows = []
    for row_idx, row in enumerate(protein_group_ids):
        is_histone = False
        protein_ids = row.strip().split(';')
        for prot_id in protein_ids:
            if prot_id in organism.histone_ids:
                is_histone = True
        if is_histone:
            histone_rows.append(row_idx)

    return histone_rows


def find_custom_proteins(protein_group_ids, custom_protein_list):
    """
    Finds those protein groups that represent the custom proteins given
    Parameters:
        - list of protein group IDs (strings of semicolon-separated IDs)
        - your list of custom proteins IDs
    Returns the row indices that represent the custom proteins
    """

    custom_protein_rows = []
    for row_idx, row in enumerate(protein_group_ids):
        is_custom_prot = False
        protein_ids = row.strip().split(';')
        for prot_id in protein_ids:
            if prot_id in custom_protein_list:
                is_custom_prot = True
        if is_custom_prot:
            custom_protein_rows.append(row_idx)

    return custom_protein_rows


def main():
    output_df, suppl_tables = process_data()
    output_df.to_csv(output_file, sep='\t', na_rep='N/A', header=True, index=False, encoding='utf-8')
    if suppl_tables:
        for idx, suppl_table in enumerate(suppl_tables):
            out_file_name, ext = os.path.splitext(output_file)
            suppl_file = "{0}_suppl_table_{1}{2}".format(out_file_name, idx, ext)
            suppl_table.to_csv(suppl_file, sep='\t', na_rep='N/A', header=True, index=True, encoding='utf-8')
    return


if __name__ == '__main__':
    setup()
    main()


####>Revision history<####
# 1.0.0 Created (VL 26/07/2019)