#!/usr/bin/env Rscript

################################################################################
# aLFQ.R         1.0.5                                                         #
# Authors: Victor Laigle (Institut Curie)                                      #
# Contact: myproms@curie.fr                                                    #
# Absolute Label-Free Quantification with R package aLFQ                       #
# Rosenberger et al., aLFQ: An R-package for estimating absolute protein       #
# quantities from label-free LC-MS/MS proteomics data. Bioinformatics. 2014    #
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

# Command line arguments and setup
cmd_args <- commandArgs(trailingOnly = TRUE)
abs_dir <- cmd_args[1]  # The directory in which (most of) the absolute quantif takes place
R_scripts_dir <- cmd_args[2]  # Needed to get FunctionLimma for normalization

source(paste0(R_scripts_dir,"/FunctionLimma.R"))

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(gtools)
library(seqinr)  # To parse fasta file
library(aLFQ)

info_file <- file.path(abs_dir, "aLFQ_info.txt")
abs_quantif_model <- file.path(abs_dir, "aLFQ_model.png")
pqi_norm_graph <- file.path(abs_dir, "pqi_norm_graph.png")
out_abs <- file.path(abs_dir, "aLFQ_abs.txt")
out_pqi<- file.path(abs_dir, "aLFQ_pqi.txt")
out_total <- file.path(abs_dir, "aLFQ_abs_total.txt")

setwd(abs_dir)

print(sessionInfo())

read_params <- function(file) {
  # Take the aLFQ_info.txt file with aLFQ parameters and read it.
  # It must be a tsv file with parameter name / parameter value pairs.
  # Some parameters are specific to aLFQ, some others are needed for this script (see below).
  # Parameters needed by aLFQ and not present in the file are set with the default value
  # (see aLFQ package documentation for more on parameters and their default values).

  # Needed arguments and parameters from info file
  info_lines <- readLines(file)
  param_names <- character(length(info_lines))
  param_values <- character(length(info_lines))

  for (i in 1:length(info_lines)) {
    info <- strsplit(info_lines[i], split = '\t')[[1]]
    param <- info[1]
    value <- info[-1]

    param_names[i] <- param
    param_values[i] <- value
  }
  params <- as.list(param_values)
  names(params) <- param_names

  # Default values for boolean parameters
  logical_params <- list(
    # Custom parameters, just in case they are not recognised
    "has_quantity" = FALSE,
    "is_mass" = FALSE,
    "is_concentration" = FALSE,
    # import function
    "averageruns" = FALSE,
    "sumruns" = FALSE,
    "openswath_superimpose_identifications" = FALSE,
    "openswath_replace_run_id" = FALSE,
    "openswath_filtertop" = FALSE,
    "openswath_removedecoys" = TRUE,
    # ProteinInference function
    "combine_precursors" = FALSE,
    "combine_peptide_sequences" = FALSE,
    "consensus_proteins" = TRUE,
    "consensus_peptides" = TRUE,
    "consensus_transitions" = TRUE
  )
  # Default values for numeric parameters
  numerical_params <- list(
    # Custom parameters, just in case they are not recognised
    "conf_level" = 0.95,
    # import function
    "mprophet_cutoff" = 0.01,
    "peptideprophet_cutoff" = 0.95,
    # ProteinInference function
    "peptide_topx" = 2,
    "transition_topx" = 3
  )
  # Default values for string parameters
  string_params <- list(
    # import function
    "abacus_column" = "ADJNSAF",
    "pepxml2csv_runsplit" = "~",
    # ProteinInference function
    "peptide_method" = "top",
    "peptide_strictness" = "strict",
    "peptide_summary" = "mean",
    "transition_strictness" = "strict",
    "transition_summary" = "sum",
    "fasta" = NA,
    "apex_model" = NA
  )

  for (param in names(logical_params)) {
    if (param %in% names(params)) {
      params[[param]] <- as.logical(params[[param]])
    } else {
      params[[param]] <- logical_params[[param]]
    }
  }
  for (param in names(numerical_params)) {
    if (param %in% names(params)) {
      params[[param]] <- as.numeric(params[[param]])
    } else {
      params[[param]] <- numerical_params[[param]]
    }
  }
  for (param in names(string_params)) {
    if (! param %in% names(params)) {
      params[[param]] <- string_params[[param]]
    }
  }

  ## Post processing of params
  if (params$model == "calibration") {
    if (is.null(params$quantities)) {  # Should never happen, need quantities for calibration
      stop("Cannot use the calibration model as required because you did not provide a file with known quantities. Make sure to provide known quantities for calibration protein or switch to proportionality model.")
    }
    params$total_protein_concentration <- 1  # Set total protein quantity to 1 with calibration model
    params$masses <- NA
  } else {  # Proportionality model
    if (is.null(params$quantities)) {  # Total quantity in each sample not provided, only percentages can be computed
      params$total_protein_concentration <- 1
      params$quantities <- NA  # Quantities are used only to calibrate the quantification
      params$masses <- NA
    } else {
      if (params$has_quantity && params$is_mass) {  # Set total quantities to 1 for aLFQ and reprocess with actual values later for masses
        params$masses <- params$quantities
        params$total_protein_concentration <- 1
        params$quantities <- NA  # Quantities are used only to calibrate the quantification
      } else {
        params$total_protein_concentration <- params$quantities
        params$quantities <- NA  # Quantities are used only to calibrate the quantification
        params$masses <- NA
      }
    }
  }
  # Coherence between parameters
  if (! params$peptide_method %in% c("iBAQ", "APEX", "NSAF")) {
    params$fasta <- string_params$fasta
  }
  if (! params$peptide_method == "APEX") {
    params$apex_model <- string_params$apex_model
  }

  return(params)
}


reformat_resultsPep <- function(resultsPep, openMS_file) {
  # Parse resultsPep.txt to extract XIC values of peptides,
  # Transform values from log back to natural intensity and
  # Build the matrix of XICs with proteins/peptides as index and states/replicates as columns.
  # This function is needed if we compute absolute quantities from our current (09/20) DDA pipeline
  # but not with DIA as long as we use OpenSwath for DIA, which is supported by aLFQ.

  keep_cols <- c('Condition', 'replicate', 'repTech', 'ProteinID', 'Peptide', 'PeptideId', 'log2Measure', 'out')
  # Computation crashes for iBAQ if ProteinID is numeric, so make it "character" explicitly. Same for "out" column.
  results <- read_tsv(resultsPep, col_types = "ccccciccdc")[ , keep_cols]

  openMS_matrix <- results %>%
    filter(is.na(out) | out == "outBalanceStates") %>%  # Exclude other outliers (usually outMixBioPeptide)
    mutate(XIC = 2 ** log2Measure, log2Measure = NULL, out = NULL) %>%
    unite("Sample", Condition, replicate, repTech, sep = '_', remove = TRUE) %>%
    separate(Peptide, c('peptide', 'charge'), sep = "_") %>%
    # Here, "peptide" is actually "seq + modifs" (e.g. PEPTIDE, but also potentially PEPTIDE&T:12).
    # Since modified peptides are considered different for quantification, it is the desired behavior
    # The modified peptides pass silently as unique peptides even with "&" and ":" symbols
    mutate(n_proteins = rep(1, nrow(.)), PeptideId = NULL) %>%
    # For now, each peptide is assigned only to one protein (head of match group).
    # TO CHANGE if we start using the information given by shared peptides, such as in the SCAMPI approach
    # Also, the same peptide is found only once per sample so PeptideId does not carry any information needed here
    spread(Sample, XIC) %>%  # "pivot_wider" used in newer versions of tidyr (1.0.0+), replaces spread
    mutate_if(is.double, round, 0) %>%
    dplyr::rename(protein = ProteinID) %>%  # Set columns names to fit the openMS format expected by aLFQ
    select(peptide, protein, n_proteins, charge, starts_with("State"))  # Same thing with columns order

  write.csv(openMS_matrix, openMS_file, quote = FALSE, row.names = FALSE)
}


aLFQ_computation <- function(params) {
  # Import primary quantification data into a unified data structure
  quantif_data <- aLFQ::import(
    ms_filenames = params$file,
    ms_filetype = params$input_type,
    concentration_filename = params$quantities,
    averageruns = params$averageruns,
    sumruns = params$sumruns,
    mprophet_cutoff = params$mprophet_cutoff,
    openswath_superimpose_identifications = params$openswath_superimpose_identifications,
    openswath_replace_run_id = params$openswath_replace_run_id,
    openswath_filtertop = params$openswath_filtertop,
    openswath_removedecoys = params$openswath_removedecoys,
    peptideprophet_cutoff = params$peptideprophet_cutoff,
    abacus_column = params$abacus_column,
    pepxml2csv_runsplit = params$pepxml2csv_runsplit
  )
  quantif_data$protein_id <- as.character(quantif_data$protein_id)  # Computation crashes for iBAQ if protein ID is numeric
  quantif_data$peptide_sequence <- gsub("&.+$", "", quantif_data$peptide_sequence, perl = TRUE)  # remove modifs in pep seq, not understood by aLFQ

  if (params$peptide_method == 'iBAQ') {
    # Parse fasta file and remove proteins not present in fasta (no sequence found in myproms database)
    # since it yields infinite iBAQ and makes every other molar percentage equal to 0.
    # Proteins without sequences in the DB might be contaminants not parsed properly at import
    sequences <- read.fasta(file = params$fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE)
    fasta_prots <- names(sequences)
    quantif_data <- quantif_data[quantif_data$protein_id %in% fasta_prots, ]
  }

  # Inference of protein quantities (actually the Protein Quantification Index values)
  prot_quantif_idx <- aLFQ::ProteinInference(
    data = quantif_data,
    peptide_method = params$peptide_method,
    peptide_topx = params$peptide_topx,
    peptide_strictness = params$peptide_strictness,
    peptide_summary = params$peptide_summary,
    transition_topx = params$transition_topx,
    transition_strictness = params$transition_strictness,
    transition_summary = params$transition_summary,
    fasta = params$fasta,
    apex_model = params$apex_model,
    combine_precursors = params$combine_precursors,
    combine_peptide_sequences = params$combine_peptide_sequences,
    consensus_proteins = params$consensus_proteins,
    consensus_peptides = params$consensus_peptides,
    consensus_transitions = params$consensus_transitions
  )

  # Check for infinite values. It can happen with iBAQ when a protein has 0 theoretical peptides
  # (according to the definition of a theoretical peptide by the iBAQ method). Not possible to check
  # it beforehand, it would require re-coding the function used in the aLFQ package.
  # If an infinite iBAQ is present, every molar % becomes 0, which is not the expected behavior
  # In that case, redo the computation after discarding the protein(s) that has an infinite PQI
  if (any(is.infinite(prot_quantif_idx$response))) {
    print("Proteins discarded because of infinite PQI :")
    print(prot_quantif_idx[is.infinite(prot_quantif_idx$response), "protein_id"])

    quantif_data <- quantif_data[quantif_data$protein_id %in% prot_quantif_idx[!is.infinite(prot_quantif_idx$response), "protein_id"], ]

    prot_quantif_idx <- aLFQ::ProteinInference(
      data = quantif_data,
      peptide_method = params$peptide_method,
      peptide_topx = params$peptide_topx,
      peptide_strictness = params$peptide_strictness,
      peptide_summary = params$peptide_summary,
      transition_topx = params$transition_topx,
      transition_strictness = params$transition_strictness,
      transition_summary = params$transition_summary,
      fasta = params$fasta,
      apex_model = params$apex_model,
      combine_precursors = params$combine_precursors,
      combine_peptide_sequences = params$combine_peptide_sequences,
      consensus_proteins = params$consensus_proteins,
      consensus_peptides = params$consensus_peptides,
      consensus_transitions = params$consensus_transitions
    )
  }

  # Format total protein concentration for each sample
  params$sample_concentration <- list()
  if (is.numeric(params$total_protein_concentration)) {
    for (sample in unique(prot_quantif_idx$run_id)) {
      params$sample_concentration[[sample]] <- params$total_protein_concentration
    }
  } else {
    concentration_data <- read.csv(params$total_protein_concentration, stringsAsFactors = FALSE)
    for (sample in concentration_data[, "run_id"]) {
      params$sample_concentration[[sample]] <- concentration_data[concentration_data$run_id == sample, "concentration"]
    }
  }

  abs_quantif <- NULL
  real_data <- NULL
  calibration_data <- tibble(
    Sample = character(),
    lm_slope = numeric(),
    lm_intercept = numeric(),
    MFE = numeric(),
    R_squared = numeric(),
    calibration_covar = numeric()
  )

  for (sample in unique(prot_quantif_idx$run_id)) {
    sample_quantif <- NULL
    sample_quantif <- aLFQ::AbsoluteQuantification(
      data = prot_quantif_idx[prot_quantif_idx$run_id == sample, ],
      total_protein_concentration = params$sample_concentration[[sample]]
    )

    if (params$model == "calibration") {
      sample_quantif <- predict(sample_quantif)
      abs_quantif <- rbind(abs_quantif, sample_quantif$prediction)
      real_data <- rbind(real_data, sample_quantif$calibration)
      calibration_data <- calibration_data %>%
        add_row(
          Sample = sample,
          lm_slope = sample_quantif$model$coefficients["response"],
          lm_intercept = sample_quantif$model$coefficients["(Intercept)"],
          MFE = sample_quantif$mfe,
          R_squared = sample_quantif$r.squared,
          calibration_covar = sample_quantif$calibration_covar
        )
    } else {
      abs_quantif <- rbind(abs_quantif, sample_quantif$estimation)
    }
  }

  calibration_data <- calibration_data %>%
    separate(Sample, c("Condition", "replicate", "repTech"), sep = "_")

  # Process normalized log output from AbsoluteQuantification function
  aLFQ_values <- NULL
  if (params$model == "proportionality") {
    aLFQ_values <- as_tibble(abs_quantif) %>%
      dplyr::left_join(as_tibble(prot_quantif_idx)[, c("run_id", "protein_id", "response")], by = c("run_id", "protein_id")) %>%
      mutate(response = round(response, 0),
             Mol_percent = signif(exp(normalized_response) * 100, 5),
             Mol_quantity = signif(exp(normalized_concentration), 5),
             normalized_response = NULL,
             normalized_concentration = NULL
      )
  } else if (params$model == "calibration") {
    aLFQ_values <- as_tibble(abs_quantif) %>%
      mutate(response = round(exp(response), 0),
             Mol_percent = signif(exp(normalized_response) * 100, 6),
             Mol_quantity = signif(exp(concentration), 6),
             normalized_response = NULL,
             concentration = NULL
      )
    real_data <- as_tibble(real_data) %>%
      separate(run_id, c("Condition", "replicate", "repTech"), sep = "_") %>%
      mutate(Prot_Quantif_Idx = exp(response),
             real_concentration = exp(concentration),
             response = NULL,
             concentration = NULL,
             normalized_response = NULL
      ) %>%
      dplyr::rename(ProteinID = protein_id)
  } else {
    stop("There was a problem in the parameters' settings, the model chosen does not exist !")
  }

  ### Compute mass values : mass% first, then mass amount (or concentration) 
  aLFQ_values <- convert_mol_to_mass(aLFQ_values, params$molecular_weights)
  
  # Proportionality model + total mass per sample provided by user:
  # - Read quantity file set aside at the beginning (total mass per sample)
  # - Compute values of Mass_quantity (with mass% and total)
  # - Recompute actual values of Mol_quantity from Mass_quantity (mol% does not change)
  if (params$model == "proportionality" && params$has_quantity && params$is_mass) {
    total_mass_data <- read_csv(params$masses, col_types = "cd")
    aLFQ_values <- aLFQ_values %>%
      dplyr::left_join(total_mass_data, by = "run_id") %>%
      group_by(run_id) %>%
      # The "concentration" column actually corresponds to the total protein mass (or total mass concentration)
      # Misnamed to avoid making a particular case earlier in the process and has no other impact
      # In all the other cases, the file is sent to an aLFQ function which requires a "concentration" column
      mutate(Mass_quantity = Mass_percent * concentration / 100, concentration = NULL) %>%
      mutate(Mol_quantity = Mass_quantity / MW) %>%
      ungroup() %>%
      select(-MW)
  } else {  # All other cases : just convert from moles to mass
    aLFQ_values <- aLFQ_values %>%
      mutate(Mass_quantity = Mol_quantity * MW) %>%
      select(-MW)
  }
  
  aLFQ_values <- aLFQ_values %>%
    separate(run_id, c("Condition", "replicate", "repTech"), sep = "_") %>%
    dplyr::rename(ProteinID = protein_id, Prot_Quantif_Idx = response) %>%
    select(Condition, replicate, repTech, everything())

  # Check if variability needs to be assessed for biological or technical replicates
  # Only one is assessed and priority is given to biological replicates
  variability <- NULL
  if (aLFQ_values %>% select(replicate) %>% n_distinct(na.rm = TRUE) > 1) {
    variability <- "biological"
  } else if (aLFQ_values %>% select(repTech) %>% n_distinct(na.rm = TRUE) > 1) {
    variability <- "technical"
  }

  # Aggregate values per condition and convert quantity units to usual ones (which will be displayed)
  mol_multiplier = 1
  mass_multiplier = 1
  if (params$is_concentration) {
    mol_multiplier = 1  # mol/L -> mmol/mL (equal, no conversion)
    mass_multiplier = 10**3  # g/L -> µg/mL (= ng/µL)
  } else {
    mol_multiplier = 10**3  # mol -> mmol
    mass_multiplier = 10**9  # g -> ng
  }
  
  # Take the geometric mean of technical replicates, then of biological replicates.
  # Geometric mean chosen over the median because there are generally only a few replicates and choosing the median
  # becomes equivalent to taking either the middle replicate or the arithmetic mean of the two middle replicates.
  # Besides, values may span several orders of magnitude, which justifies the geometric mean.
  abs_values <- aLFQ_values %>%
    dplyr::select(Condition, replicate, repTech, ProteinID, Mol_quantity, Mass_quantity) %>%
    group_by(ProteinID, Condition, replicate) %>%
    dplyr::summarise(
      M_geo_sd = if (!is.null(variability) && variability == "technical") geo_sd(Mol_quantity) else NA,  # same mol/mass
      M_geo_cv = if (!is.null(variability) && variability == "technical") geo_cv(Mol_quantity) else NA,  # same mol/mass
      Mol_CIs  = if (!is.null(variability) && variability == "technical") list(geo_ci(Mol_quantity, conf.level = params$conf_level)) else NA,
      Mass_CIs = if (!is.null(variability) && variability == "technical") list(geo_ci(Mass_quantity, conf.level = params$conf_level)) else NA,
      Mol_quantity  = geo_mean(Mol_quantity),  # Keep means at the end, the order matters
      Mass_quantity = geo_mean(Mass_quantity)  # Keep means at the end, the order matters
    ) %>%
    ungroup()
  
  abs_values <- abs_values %>%
    group_by(ProteinID, Condition) %>%
    dplyr::summarise(  # If variability is technical, only 1 bio rep, so max/first is evaluated on a single value (could have been anything)
      M_geo_sd = if (is.null(variability)) NA else if (variability == "technical") max(M_geo_sd) else geo_sd(Mol_quantity),  # same mol/mass
      M_geo_cv = if (is.null(variability)) NA else if (variability == "technical") max(M_geo_cv) else geo_cv(Mol_quantity),  # same mol/mass
      Mol_CIs  = if (is.null(variability)) list(tibble(ci_inf = NA, ci_sup = NA)) else if (variability == "technical") list(first(Mol_CIs)) else list(geo_ci(Mol_quantity, conf.level = params$conf_level)),
      Mass_CIs = if (is.null(variability)) list(tibble(ci_inf = NA, ci_sup = NA)) else if (variability == "technical") list(first(Mass_CIs)) else list(geo_ci(Mass_quantity, conf.level = params$conf_level)),
      Mol_quantity  = geo_mean(Mol_quantity),  # Keep means at the end, the order matters
      Mass_quantity = geo_mean(Mass_quantity)  # Keep means at the end, the order matters
    ) %>%
    unnest(Mol_CIs, Mass_CIs) %>%
    ungroup()

  # Recompute percentages because we cannot directly take the mean on percentages (for varying total quantities).
  # In the calibration case, the accuracy of the percentages really depends on the accuracy of the model.
  # We also report 95% CIs for the percentages means, but only in the proportionality case. Proper error estimation
  # is not obvious at all because of the ratios and sums. However, the sums are constant in the proportionality case,
  # which means that errors on individual proteins balance one another. This assumption of constant sums in the
  # proportionality case is always verified, whether the total quantity of proteins is given by the user (sum known for
  # either moles or mass) or not (sums to 1). In the calibration case, the CIs are not computed (yet) because such an
  # assumption cannot be verified.
  if (params$model == "proportionality") {
    abs_values <- abs_values %>%
      rename(Mol_CI_inf = ci_inf, Mol_CI_sup = ci_sup, Mass_CI_inf = ci_inf1, Mass_CI_sup = ci_sup1) %>%
      group_by(Condition) %>%
      mutate(
        Mol_percent         = signif(Mol_quantity / sum(Mol_quantity) * 100, 5),
        Mol_percent_CI_inf  = signif(Mol_CI_inf / sum(Mol_quantity) * 100, 5),
        Mol_percent_CI_sup  = signif(pmin(100, Mol_CI_sup / sum(Mol_quantity) * 100), 5),  # Cap CI_sup to 100%
        Mass_percent        = signif(Mass_quantity / sum(Mass_quantity) * 100, 5),
        Mass_percent_CI_inf = signif(Mass_CI_inf / sum(Mass_quantity) * 100, 5),
        Mass_percent_CI_sup = signif(pmin(100, Mass_CI_sup / sum(Mass_quantity) * 100), 5)  # Cap CI_sup to 100%
      ) %>%
      ungroup()
  } else {
    abs_values <- abs_values %>%
      rename(Mol_CI_inf = ci_inf, Mol_CI_sup = ci_sup, Mass_CI_inf = ci_inf1, Mass_CI_sup = ci_sup1) %>%
      group_by(Condition) %>%
      mutate(
        Mol_percent         = signif(Mol_quantity / sum(Mol_quantity) * 100, 5),
        Mol_percent_CI_inf  = NA,
        Mol_percent_CI_sup  = NA,
        Mass_percent        = signif(Mass_quantity / sum(Mass_quantity) * 100, 5),
        Mass_percent_CI_inf = NA,
        Mass_percent_CI_sup = NA
      ) %>%
      ungroup()
  }

  abs_values <- abs_values %>%
    mutate(  # convert quantities and CIs to relevant units 
      Mol_quantity  = signif(Mol_quantity * mol_multiplier, 5),
      Mol_CI_inf    = signif(Mol_CI_inf * mol_multiplier, 5),
      Mol_CI_sup    = signif(Mol_CI_sup * mol_multiplier, 5),
      Mass_quantity = signif(Mass_quantity * mass_multiplier, 5),
      Mass_CI_inf   = signif(Mass_CI_inf * mass_multiplier, 5),
      Mass_CI_sup   = signif(Mass_CI_sup * mass_multiplier, 5)
    ) %>%
    select(  # Reorder columns
      Condition, ProteinID,
      Mol_quantity, Mol_CI_inf, Mol_CI_sup,
      Mol_percent, Mol_percent_CI_inf, Mol_percent_CI_sup,
      Mass_quantity, Mass_CI_inf, Mass_CI_sup,
      Mass_percent, Mass_percent_CI_inf, Mass_percent_CI_sup,
      everything()
    )

  return(list(abs = abs_values, pqi = aLFQ_values, calib = calibration_data, real = real_data))
}


plot_model <- function(aLFQ_matrix, model_file, calib_data = NULL, real_data = NULL) {
  # Plot the transformation of PQI to absolute quantities for each sample.
  # If anchor proteins were used for calibration it adds the real values for these proteins,
  # and the linear model, R squared and Mean Fold Error computed by aLFQ.
  quantif_data <- aLFQ_matrix %>%
    mutate(Condition = factor(Condition, levels = mixedsort(unique(Condition))),
           replicate = factor(replicate, levels = mixedsort(unique(replicate))),
           repTech = factor(repTech, levels = mixedsort(unique(repTech)))
          )
  n_samples <- n_distinct(aLFQ_matrix$Condition, aLFQ_matrix$replicate, aLFQ_matrix$repTech)
  fig_size <- ifelse(ceiling(sqrt(n_samples)) * 2 > 15, 15, ceiling(sqrt(n_samples)) * 2)

  calib_plot <- ggplot(quantif_data) +
    geom_point(aes(x = log(Prot_Quantif_Idx), y = log(Mol_quantity)), colour = "darkblue", alpha = 0.4, size = 0.01)
  if (!is.null(calib_data)) {
    model_data <- calib_data %>%
      mutate(Condition = factor(Condition, levels = mixedsort(unique(Condition))),
             replicate = factor(replicate, levels = mixedsort(unique(replicate))),
             repTech = factor(repTech, levels = mixedsort(unique(repTech)))
            )
    annot_lm <- tibble(Condition = model_data$Condition,
                         replicate = model_data$replicate,
                         repTech = model_data$repTech,
                         text = paste0("y = ", round(model_data$lm_slope, 3), "x+", round(model_data$lm_intercept, 2), "\n R2 = ", round(model_data$R_squared, 3))
                         )
    annot_error <- tibble(Condition = model_data$Condition,
                         replicate = model_data$replicate,
                         repTech = model_data$repTech,
                         text = paste0("MFE = ", round(model_data$MFE, 2))
                         )
    calib_plot <- calib_plot +
      geom_abline(data = model_data, aes(slope = lm_slope, intercept = lm_intercept), color = "red", size = 0.1) +
      geom_text(data = annot_lm, aes(x = -Inf, y = Inf, label = text, hjust = -0.1, vjust = 1.1), size = fig_size / 3) +
      geom_text(data = annot_error, aes(x = Inf, y = -Inf, label = text, hjust = 1.1, vjust = -0.5), size = fig_size / 3)
  }
  if (!is.null(real_data)) {
    real <- real_data %>%
      mutate(Condition = factor(Condition, levels = mixedsort(unique(Condition))),
             replicate = factor(replicate, levels = mixedsort(unique(replicate))),
             repTech = factor(repTech, levels = mixedsort(unique(repTech)))
            )
    calib_plot <- calib_plot +
      geom_point(data = real, aes(x = log(Prot_Quantif_Idx), y = log(real_concentration)), colour = "darkred", size = 0.5)
  }
  calib_plot <- calib_plot +
    theme(strip.text = element_text(size = fig_size / 1.1)) +
    facet_wrap(vars(Condition, replicate, repTech), scale = "fixed", labeller = label_wrap_gen(multi_line = FALSE)) +
    ylab("Log(Absolute Quantity)") +
    xlab("Log(Protein Quantification Index)")

  ggsave(model_file, plot = calib_plot, width = fig_size, height = fig_size)
}


normalize_PQI <- function(pqi_data, normalization_type) {
  # Function to normalize the Protein Quantification Indices across the different samples
  # to make the PQI values directly comparable. Basically, we use a robust scaler on each sample
  # (substract median and divide by MAD) and we re-scale the values to the approximate original range
  # with the mean of the medians and the geometric mean of the MADs

  # Format dataframe as expected by normalization function from FunctionLimma.R
  # The (little) trick is to use the expected "experiment" column (not used currently) in an artificial way
  # to normalize either between states or wihin states
  to_norm <- NULL
  if (normalization_type == "global") {
    to_norm <- pqi_data %>%
      dplyr::rename(sample = Condition, M = Prot_Quantif_Idx) %>%
      mutate(experiment = "A", normProtein = 1, peptideCom = 1)  # Same experiment name to normalize across all states
  } else if (normalization_type == "per_state") {
    to_norm <- pqi_data %>%
      dplyr::rename(sample = Condition, M = Prot_Quantif_Idx) %>%
      mutate(experiment = sample, normProtein = 1, peptideCom = 1)  # Different experiment names to normalize within states
  }

  # Call to normalization function from FunctionLimma.R with median centering and MAD scaling
  centering <- median
  scaling <- mad
  normalized_data <- NULL
  if (normalization_type == "none") {
    normalized_data <- pqi_data
  } else {
    X <- .normWithinExp(to_norm, centering, scaling, "PEP_INTENSITY")
    # Back to original shape of data
    normalized_data <- X$data %>%
      dplyr::rename(Condition = sample, Prot_Quantif_Idx = M) %>%
      mutate(Prot_Quantif_Idx = round(Prot_Quantif_Idx, 0), experiment = NULL, normProtein = NULL, peptideCom = NULL)
  }

  return(normalized_data)
}


plot_PQI_normalization <- function(old_pqi, new_pqi, norm_file, normalization_type) {
  plot_title <- NULL
  pqi_matrices <- NULL
  if (normalization_type == "global") {
    plot_title <- "PQI before and after robust normalization (median and scale) on all samples"
    pqi_matrices <- list("Before normalization" = old_pqi, "After normalization" = new_pqi)
  } else if (normalization_type == "per_state") {
    plot_title <- "PQI before and after robust normalization (median and scale) within conditions"
    pqi_matrices <- list("Before normalization" = old_pqi, "After normalization" = new_pqi)
  } else {
    plot_title <- "PQI without normalization"
    pqi_matrices <- list("Without normalization" = old_pqi)
  }

  pqi_data <- bind_rows(pqi_matrices, .id = "Normalization") %>%
    mutate(Sample = paste(Condition, replicate, repTech, sep = "_"),
           Normalization = factor(Normalization, levels = sort(unique(Normalization), decreasing = TRUE))  # decreasing to ensure that "before" is actually displayed before
    ) %>%
    mutate(Sample = factor(Sample, levels = mixedsort(unique(Sample))))

  pqi_norm_plot <- ggplot(pqi_data) +
    aes(x = Sample, y = log(Prot_Quantif_Idx), color = Condition) +
    geom_boxplot() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_wrap(vars(Normalization)) +
    labs(title = plot_title)

  ggsave(norm_file, plot = pqi_norm_plot, width = 12, height = 6)
}


convert_mol_to_mass <- function(data, mw_file) {
  # Convert molar percentages to mass percentages for the whole set of proteins
  # data must contain at least 3 columns : protein ID, run_id, Mol_percent
  # mw_file is a tsv file containing the myproms ID_PROTEIN of each protein and its respective MW
  
  mol_weights <- read_tsv(mw_file, col_types = "cd")  # Load molecular weights of proteins
  
  # Convert mol% to mass%
  new_data <- data %>%
    dplyr::inner_join(mol_weights, by = c("protein_id" = "ID_PROTEIN")) %>%  # discard protein if no MW
    group_by(run_id) %>%
    mutate(Mass_percent = 100 * Mol_percent * MW / sum(Mol_percent * MW)) %>%
    ungroup()

  return(new_data)
}


geo_mean <- function(x) {
  return(exp(mean(log(x), na.rm = TRUE)))
}

geo_sd <- function(x) {
  return(exp(sd(log(x), na.rm = TRUE)))
}

geo_cv <- function(x) {
  return(sqrt(exp(sd(log(x), na.rm = TRUE)^2) - 1) * 100)
}

geo_ci <- function(x, conf.level = 0.95, alternative = "two_tailed") {
  if (conf.level <= 0 || conf.level >= 1) {
    stop(paste0(c("The confidence level must be between 0 and 1, received ", conf.level)))
  }
  if (length(x) < 2) {
    return(tibble(ci_inf = NA, ci_sup = NA))
  }
  m <- geo_mean(x)
  s <- geo_sd(x)
  n <- length(which(!is.na(x)))
  df <- n-1
  if(alternative == "greater") {
    quantile <- qt(conf.level, df = df)
    ci_inf <- min(m / (quantile * s/sqrt(n)), m * (quantile * s/sqrt(n)), na.rm=TRUE)
    ci_sup <- +Inf
  } else if (alternative == "less") {
    quantile <- qt(conf.level, df = df)
    ci_inf <- 0
    ci_sup <- max(m / (quantile * s/sqrt(n)), m * (quantile * s/sqrt(n)), na.rm=TRUE)
  } else{
    quantile <- qt((1 - (1-conf.level)/2), df = df)
    ci_inf <- min(m / (quantile * s/sqrt(n)), m * (quantile * s/sqrt(n)), na.rm=TRUE)
    ci_sup <- max(m / (quantile * s/sqrt(n)), m * (quantile * s/sqrt(n)), na.rm=TRUE)
  }
  return(tibble(ci_inf = ci_inf, ci_sup = ci_sup))
}


################################################################################
##### Main computation #########################################################
################################################################################

if (!file.exists(info_file)) {
  stop(paste("Information file", info_file, "does not exist, make sure that the right path was provided and that it has been written.", sep = " "))
}
params <- read_params(info_file)

# Test if absolute quantif for DDA or DIA
# - True = DDA = convert to openMS format
# - False = DIA = keep openswath input as is
if ("resultsPep" %in% names(params)) {
  openMS_file <- paste0(abs_dir, "/openMSlfq.txt")
  reformat_resultsPep(params$resultsPep, openMS_file)
  if (params$input_type == "openmslfq") {  # This should always be the case if resultsPep is present
    params$file <- openMS_file
  }
}

aLFQ_results <- aLFQ_computation(params)

if (params$model == "calibration") {
  plot_model(aLFQ_results$pqi, abs_quantif_model, aLFQ_results$calib, aLFQ_results$real)
} else {
  plot_model(aLFQ_results$pqi, abs_quantif_model)
}

normalized_pqi <- normalize_PQI(aLFQ_results$pqi, params$pqi_normalization)
plot_PQI_normalization(aLFQ_results$pqi, normalized_pqi, pqi_norm_graph, params$pqi_normalization)

# Keep quantity only if total or anchor proteins quantity was provided and
# Rename quantity according to the provided units
abs_total <- NULL
if (params$has_quantity) {
  if (params$is_concentration) {
    aLFQ_results$abs <- aLFQ_results$abs %>%
      rename(Mol_concentration = Mol_quantity, Mass_concentration = Mass_quantity)
    abs_total <- aLFQ_results$abs %>%
      group_by(Condition) %>%
      summarise(Mol_concentration = sum(Mol_concentration), Mass_concentration = sum(Mass_concentration))
  } else {
    aLFQ_results$abs <- aLFQ_results$abs %>%
      rename(Mol_amount = Mol_quantity, Mass_amount = Mass_quantity)
    abs_total <- aLFQ_results$abs %>%
      group_by(Condition) %>%
      summarise(Mol_amount = sum(Mol_amount), Mass_amount = sum(Mass_amount))
  }
} else {  # Remove variables that correspond to actual quantities, keep only percentages
  aLFQ_results$abs <- aLFQ_results$abs %>%
    mutate(Mol_quantity = NULL, Mol_CI_inf = NULL, Mol_CI_sup = NULL,
           Mass_quantity = NULL, Mass_CI_inf = NULL, Mass_CI_sup = NULL
    )
}

write.table(aLFQ_results$abs, out_abs, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(normalized_pqi, out_pqi, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
if (params$has_quantity) {
  write.table(abs_total, out_total, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

####>Revision history<####
# 1.0.5 [ENHANCEMENT] Add computation of geometric SD, CV and 95% CIs for variability and error estimation (VL 07/12/20)
# 1.0.4 [ENHANCEMENT] Add computation of total amount per state + [BUGFIX] Handle modified peptides (VL 01/12/20)
# 1.0.3 [BUGFIX] Change tibble func to as_tibble in some cases for compatibility with downgraded versions (VL 27/11/20)
# 1.0.2 [FEATURE] Add handling of mass quantities (VL 13/10/20)
# 1.0.1 [ENHANCEMENT] Add handling of known protein quantities, normalization of PQI, plot of the normalization and plot of the correspondance between PQI and absolute quantities (VL 04/09/20)
# 1.0.0 Created (VL 19/08/2020)
