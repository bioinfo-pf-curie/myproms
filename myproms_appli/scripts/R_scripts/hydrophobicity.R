#!/usr/bin/env Rscript

################################################################################
# hydrophobicity.R         1.0.0                                               #
# Authors: Victor Laigle (Institut Curie)                                      #
# Contact: myproms@curie.fr                                                    #
# Computation of peptides hydrophobicity based on SSRCalc algorithm            #
# Krokhin et al., An improved model for prediction of retention times of       #
# tryptic peptides in ion pair reversed-phase HPLC: its application to protein #
# peptide mapping by off-line HPLC-MALDI MS, Mol. Cell. Proteomics, 2004.      #
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

library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ggplot2)
library(hexbin)
library(scales)  # For hue_pal


##########################################
### Version 1 of the SSRCalc algorithm ###
##########################################

# Step 1: Peptides retention coefficients for C18 100 Å column and Formic Acid
# Dwivedi et al., Practical implementation of 2D HPLC scheme with accurate peptide retention
# prediction in both dimensions for high-throughput bottom-up proteomics, Anal. Chem., 2008.

# Step 2: Equation to get peptides retention coefficients for N terminal residues
# Step 3: Equation and coefficients to compute Hydrophobicity value (H)
# Step 4: Correction coefficient (KL) to account for peptide length (N)
# Step 5: Correction factor based on computed hydrophobicity (to correct for high values)
# Krokhin et al., An improved model for prediction of retention times of tryptic peptides
# in ion pair reversed-phase HPLC: its application to protein peptide mapping by off-line
# HPLC-MALDI MS, Mol. Cell. Proteomics, 2004.

# Step 6: Equation to switch from computed hydrophobicity to HI with C18 100 Å column and Formic Acid
# Krokhin and Spicer, Peptide retention standards and hydrophobicity indexes in reversed-phase 
# high-performance liquid chromatography of peptides, Anal. Chem., 2009.


# Step 1: retention coefficients
rt_coeff <- tribble(
  ~Amino_Acid, ~Rc,
  "W", 13.67,
  "F", 11.92,
  "L", 9.89,
  "I", 9.06,
  "M", 6.96,
  "V", 5.72,
  "Y", 5.97,
  "C", 0.7,
  "P", 1.98,
  "A", 1.63,
  "E", 1.75,
  "T", 1.37,
  "D", 0.95,
  "Q", 0.2,
  "S", 0.27,
  "G", -0.07,
  "R", -3.55,
  "N", -0.59,
  "H", -5.05,
  "K", -5.08
)

# Step 2: N-terminal retention coefficients
rt_coeff <- rt_coeff %>% 
  mutate(RcNt = sum(Rc)/20 - Rc) %>% 
  column_to_rownames("Amino_Acid")

compute_HI_v1 <- Vectorize(function(seq) {

  seq_list <- str_split(seq, "")[[1]]
  N <- length(seq_list)

  # Step 3
  rc_sum <- sum(sapply(seq_list, function(x) rt_coeff[x, "Rc"]))
  H <- rc_sum + 0.42 * rt_coeff[seq_list[1], "RcNt"] + 0.22 * rt_coeff[seq_list[2], "RcNt"] + 0.05 * rt_coeff[seq_list[3], "RcNt"]
  
  # Step 4
  KL <- case_when(
    N < 10 ~ 1 - 0.027 * (10 - N),
    N > 20 ~ 1 - 0.014 * (N - 20),
    TRUE ~ 1
  )  
  H <- H * KL
  
  # Step 5
  H_final <- case_when(
    H < 38 ~ H,
    H >= 38 ~ H - 0.3 * (H - 38)
  )

  # Step 6
  HI <- 0.4954 * H_final - 2.6687
  
  return(HI)
})


######################################################
### Version 2/3 (partial) of the SSRCalc algorithm ###
######################################################

# Step 1: Peptides retention coefficients for C18 100 Å column (but based on TFA)
# Krokhin, Sequence-specific retention calculator. Algorithm for peptide retention prediction 
# in ion-pair RP-HPLC: application to 300- and 100-A pore size C18 sorbents, Anal. Chem., 2006. 

# Step 2: Equation and coefficients to compute Hydrophobicity value (H)
# Guessed

# Step 3: Correction coefficient (KL) to account for peptide length (N)
# Step 4: Correction factor based on computed hydrophobicity (to correct for high values)
# Krokhin, Sequence-specific retention calculator. Algorithm for peptide retention prediction 
# in ion-pair RP-HPLC: application to 300- and 100-A pore size C18 sorbents, Anal. Chem., 2006. 

# Step 5: Equation to switch from computed hydrophobicity to HI with C18 100 Å column and Formic Acid
# Krokhin and Spicer, Peptide retention standards and hydrophobicity indexes in reversed-phase 
# high-performance liquid chromatography of peptides, Anal. Chem., 2009.

# Step 1 and 2
# rt_coeff_v2 <- tribble(  # Rc for peptides with N≥9, Rcs for peptides with N<9
#   ~Amino_Acid, ~Rc, ~Rc1, ~Rc2, ~RcN, ~`Rc(N-1)`, ~Rcs, ~Rc1s, ~Rc2s, ~RcNs, ~`Rc(N-1)s`,
#   "W", 13.35, 11.5, 11.8, 13.35, 13, 13.9, 11.8, 13, 13.9, 12.9,
#   "F", 11.67, 7.6, 9.7, 11.67, 11.5, 11.3, 8.4, 10, 11.3, 10.85,
#   "L", 9.4, 5.57, 7.4, 9.4, 9.3, 8.7, 5.5, 7.7, 8.7, 8.5,
#   "I", 7.96, 4.95, 6.3, 7.96, 6.6, 7.25, 4.5, 6.5, 7.25, 7.2,
#   "M", 6.27, 5.2, 5.7, 6.27, 5.8, 6.25, 4.2, 5.7, 6.25, 5.6,
#   "V", 4.68, 2.1, 3.4, 4.68, 3.9, 4.4, 2.1, 3, 4.4, 4.4,
#   "Y", 5.35, 4.3, 5.1, 5.35, 5, 5.7, 5, 5.4, 5.7, 5.3,
#   "C", 0.1, 0.4, 0.2, 0.1, -0.4, 0.6, 0.6, 1, 0.6, -0.5,
#   "P", 1.85, 1.7, 1.75, 1.85, 1.2, 2.5, 1.7, 2.1, 2.5, 1.9,
#   "A", 1.02, -0.35, 0.35, 1.02, -0.2, 0.5, -0.05, 0.1, 0.5, -0.3,
#   "E", 1, 1, -0.2, 1, -0.1, 0.7, 0.45, 0.5, 0, 0.25,
#   "T", 0.64, 0.95, 0.6, 0.64, -0.1, 0.4, 0.3, 0.4, 0.4, -0.5,
#   "D", 0.15, 0.9, 0.6, 0.15, -0.4, 0.6, 0.3, 0.2, 0.6, -0.5,
#   "Q", -0.6, -0.5, -0.2, -0.6, -1.1, -0.4, -0.2, -0.7, -0.4, -1.3,
#   "S", -0.14, 1.1, -0.1, -0.14, -1, -0.4, 0.2, -0.3, -0.4, -1.2,
#   "G", -0.35, 0.15, 0.15, -0.35, -0.4, 0, 0.15, 0.2, 0, -0.7,
#   "R", -2.55, -1.4, -1.5, -1.1, -1.3, -1, 0.4, -1, -1.1, -1.9,
#   "N", -0.95, 1.2, -0.1, -0.95, -1.3, -0.65, 0.4, -0.05, -0.65, -1.2,
#   "H", -3, -1.4, -1, -3, -1.9, -1.3, -1.3, -1.1, -1.3, -1.7,
#   "K", -3.4, -1.85, -2.3, -2.1, -2.1, -1.75, -1.5, -1.75, -2.3, -2.5
# ) %>% 
#   column_to_rownames("Amino_Acid")

# compute_HI_v2 <- Vectorize(function(seq) {
#   seq_list <- str_split(seq, "")[[1]]
#   N <- length(seq_list)
  
#   # Step 2
#   H <- NULL
#   if (N >= 9) {
#     rc_sum <- sum(sapply(seq_list[c(3:(N-2))], function(x) rt_coeff_v2[x, "Rc"]))
#     H <-  rt_coeff_v2[seq_list[1], "Rc1"] + 
#       rt_coeff_v2[seq_list[2], "Rc2"] + 
#       rc_sum + 
#       rt_coeff_v2[seq_list[N-1], "Rc(N-1)"] + 
#       rt_coeff_v2[seq_list[N], "RcN"]
#   } else {
#     rc_sum <- 0
#     if (N > 4) {
#       rc_sum <- sum(sapply(seq_list[c(3:(N-2))], function(x) rt_coeff_v2[x, "Rcs"]))
#     }
#     H <-  rt_coeff_v2[seq_list[1], "Rc1s"] + 
#       rt_coeff_v2[seq_list[2], "Rc2s"] + 
#       rc_sum + 
#       rt_coeff_v2[seq_list[N-1], "Rc(N-1)s"] + 
#       rt_coeff_v2[seq_list[N], "RcNs"]
#   }
  
#   # Step 3
#   KL <- case_when(
#     N < 8 ~ 1 - 0.055 * (8 - N),
#     N > 20 ~ 1 / (1 + 0.027 * (N - 20)),
#     TRUE ~ 1
#   )  
#   H <- H * KL
  
#   # Step 4
#   H_final <- NULL
#   if (H <= 20) {
#     H_final <- H
#   } else if (20 < H & H <= 30) {
#     H_final <- H - 0.27 * (H - 18)
#   } else if (30 < H & H <= 40) {
#     H_final <- H - 0.33 * (H - 18)
#   } else if (40 < H & H <= 50) {
#     H_final <- H - 0.38 * (H - 18)
#   } else if (H > 50) {
#     H_final <- H - 0.447 * (H - 18)
#   }
  
#   # Step 5
#   HI <- 0.4954 * H_final - 2.6687
  
#   return(HI)
# })

#######################
### Other functions ###
#######################

read_params <- function(file) {
  # Take the hydro_info.txt file with parameters for hydrophobicity computation and read it.
  # It must be a tsv file with parameter name / parameter value pairs.

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

  if (!("plot_as" %in% names(params))) {
    params[["plot_as"]] <- "png"
  }
  if (("rt_offset" %in% names(params))) {
    params[["rt_offset"]] <- as.numeric(params[["rt_offset"]])
  } else {
    params[["rt_offset"]] <- 0
  }

  return(params)
}


compute_hydro <- function(data_file, compute_HI = compute_HI_v1) {
  full_data <- read_tsv(data_file, col_types = "dcdcdcddd")  # AnaID, AnaName, SampID, SampName, PepID, Seq, RT, Charge, Status
  pep_data <- full_data %>%
    select(Sequence) %>%
    distinct(Sequence) %>%  # Compute HI only once for each sequence
    mutate(HI = compute_HI(Sequence))

  full_data <- full_data %>%
    left_join(pep_data, by = "Sequence")

  return(full_data)
}


compute_linear_models <- function(df, rt_offset = 0) {
  lin_mod <- df %>%
    filter(RT > rt_offset) %>%
    group_by(Analysis_ID, Analysis, Sample_ID, Sample) %>%
    summarise(slope = round(lm(HI ~ RT)$coefficients[2], 3),
              intercept = round(lm(HI ~ RT)$coefficients[1], 3),
              R = round(cor(HI, RT), 4),
              R2 = round(cor(HI, RT) ** 2, 4)
    ) %>%
    mutate(label = ifelse(intercept >= 0, 
                          paste0("y = ", slope, "x + ", intercept, "\nR2 = ", R2), 
                          paste0("y = ", slope, "x - ", -intercept, "\nR2 = ", R2)
                          )
    ) %>%
    ungroup()

  return(lin_mod)
}


compute_percentile_time_shift <- function(df, lm_df) {
  # For each peptide, compute the error between observed RT and predicted RT (with linear model and HI).
  # Then compute the percentile at which each peptide is (per Analysis).

  percentile_df <- df %>% 
    left_join(lm_df, by = c("Analysis_ID", "Analysis", "Sample_ID", "Sample")) %>%
    mutate(predicted_RT = round((HI - intercept) / slope, 6), time_shift = round((RT - (HI - intercept) / slope), 6)) %>%
    group_by(Analysis_ID) %>%
    mutate(percentile = round(percent_rank(abs(time_shift)), 6)) %>%
    ungroup() %>%
    select(-slope, -intercept, -R, -R2, -label)
  
  percentile_df$percentile <- format(percentile_df$percentile, scientific = FALSE)

  return(percentile_df)  
}


determine_row_col_nb <- function(n_samples, max_ana_per_samp) {
  # Try to adapt number of rows/columns of graphs to number of analyses and samples
  nrows <- NULL
  ncols <- NULL
  if (max_ana_per_samp > 1) {
    if (max_ana_per_samp <= 3) {
      if (n_samples < 5) {
        nrows <- n_samples
      } else if (n_samples < 10) {
        ncols <- max_ana_per_samp * 2
      } else if (n_samples < 30) {
        ncols <- max_ana_per_samp * 3
      } else {
        ncols <- round(sqrt(max_ana_per_samp * n_samples))
      }
    } else {
      if (n_samples < 5) {
        nrows <- n_samples
      } else {
        ncols <- round(sqrt(max_ana_per_samp * n_samples))
      }
    }          
  } else if (n_samples <= 3) {
    nrows <- 1
  } else if (n_samples <= 8) {
    nrows <- 2
  }
  return(list(nrows = nrows, ncols = ncols))
}


determine_fig_width_height <- function(nrows, ncols, n_analyses, n_samples, max_ana_per_samp) {
  fig_width <- 10
  fig_height <- 10
  if (n_analyses > 1) {
    if (n_samples > 1) {
      if (is.null(ncols) && is.null(nrows)) {  # 1 analysis per sample
        fig_width <- min(5 * ceiling(sqrt(max_ana_per_samp * n_samples)), 30)
        fig_height <- min(5 * n_samples + 5, 30) + 3
      } else if (!is.null(nrows)) {
        fig_width <- min(5 * ceiling(sqrt(max_ana_per_samp * n_samples)), 30)
        fig_height <- min(5 * nrows, 30) + 3
      } else {
        fig_width <- min(5 * ncols, 30)
        fig_height <- min(5 * ncols, 30) + 3
      }
    } else {
      fig_width <- min(5 * n_analyses, 30)
      fig_height <- min(5 * n_analyses, 30)
    }
  }

  return(list(fig_width = fig_width, fig_height = fig_height))
}


plot_HI_RT_shifts <- function(df, lm_df, plot_file, plot_as = "png", rt_offset = 0) {
  size_ratio <- 1 / 0.352777778  # Ratio of text size between geom_text and other texts (pt <-> in or similar)
  text_size <- 8
  
  df <- df %>%
    mutate(graph_label = paste(Sample, Analysis, sep = ":"),
           samp_color = hue_pal()(n_distinct(df$Sample))[as.numeric(as.factor(Sample))]
    )
  lm_df <- lm_df %>%
    mutate(graph_label = paste(Sample, Analysis, sep = ":"))
  
  strip_sample <- function(strg) {
    return(sub(".+:", "", strg))
  }
  
  # Function to change color of facet label background, to distinguish between samples
  change_facet_label_color <- function(table_grob, ana_colors) {
    label <- table_grob[["grobs"]][[1]][["children"]][[2]][["children"]][[1]][["label"]]
    if (!is.null(label)) {  # Empty facets do not have labels
      table_grob[["grobs"]][[1]][["children"]][[1]]$gp$fill <- ana_colors[ana_colors$Analysis == label, "samp_color"] %>% pull()
    }
    return(table_grob)
  }
  
  # Function that draws 1 plot, whether it is experiment, sample or analysis scope
  draw_graph <- function(df, graph_type, lm_df, plot_file, scope) {
    n_analyses <- n_distinct(df$Analysis_ID)
    n_samples <- n_distinct(df$Sample_ID)
    max_ana_per_samp <- df %>% 
      select(Analysis_ID, Sample_ID) %>% 
      group_by(Sample_ID) %>% 
      summarise(Analysis_ID = n_distinct(Analysis_ID), .groups = "drop_last") %>% 
      summarise(max = max(Analysis_ID)) %>% 
      pull()
    fig_width <- 10
    fig_height <- 10
    pep_nb <- length(df$RT)
    min_x_value <- min(df$RT)
    max_y_value <- max(df$HI)
    max_RT <- max(df$RT)
    shift_1pthd <- quantile(df$time_shift, probs = seq(0, 1, by= 0.001))["0.1%"]
    shift_999pthd <- quantile(df$time_shift, probs = seq(0, 1, by= 0.001))["99.9%"]

    # Draw main graph (points / density or distribution)
    p <- ggplot(df)
    if (graph_type == "shift_distrib") {
      p <- p +
        aes(x = time_shift) +
        geom_histogram(aes(fill = Sample), binwidth = max_RT / 500) +
        xlab("RT shift betwwen obs. and pred. (min)") +
        xlim(shift_1pthd, shift_999pthd)
    } else if (graph_type == "shift_vs_RT") {
      if (pep_nb <= 100 * n_analyses) {  # Density graphs are not that good on small datasets
        p <- p +
          aes(x = RT, y = time_shift) +
          geom_hline(yintercept = 0, size = 0.3, color = "grey40", linetype = "dashed") +
          geom_point(aes(color = Sample), size = 0.5) +
          xlab("Observed RT (min)") +
          ylab("Shift between obs. and pred. RT (min)")
      } else {
        p <- p +
          aes(x = RT, y = time_shift) +
          geom_hex(bins = 60) +  # Bins as large in minutes as the run was in hours to adapt size
          scale_fill_viridis_c(option = "plasma") +
          geom_hline(yintercept = 0, size = 0.1, color = "white", linetype = "dashed") +
          guides(fill = guide_colorbar(title.vjust = .8)) +
          xlab("Observed RT (min)") +
          ylab("Shift between obs. and pred. RT (min)")
      }
    } else {  # graph_type == "HI_vs_RT"
      if (pep_nb <= 100 * n_analyses) {  # Density graphs are not that good on small datasets
        p <- p +
          aes(x = RT, y = HI) +
          geom_abline(data = lm_df, aes(intercept = intercept, slope = slope), size = 0.5, color = "grey40", linetype = "dashed") +
          geom_vline(xintercept = rt_offset, size = 0.1, color = "black", linetype = "solid") +
          geom_point(aes(color = Sample), size = 0.5) +
          geom_text(data = lm_df, aes(x = min_x_value + 5, y = 0.9 * max_y_value, label = label), hjust = 0, fontface = "plain", size = text_size / size_ratio) +
          guides(fill = guide_colorbar(title.vjust = .8)) +
          xlab("RT (min)") +
          ylab("HI (%ACN)")
      } else {
        p <- p +
          aes(x = RT, y = HI) +
          geom_hex(bins = 60) +  # Bins as large in minutes as the run was in hours to adapt size
          scale_fill_viridis_c(option = "plasma") +
          geom_abline(data = lm_df, aes(intercept = intercept, slope = slope), size = 0.3, color = "white", linetype = "dashed") +
          geom_vline(xintercept = rt_offset, size = 0.1, color = "black", linetype = "solid") +
          geom_text(data = lm_df, aes(x = min_x_value + 5, y = 0.9 * max_y_value, label = label), hjust = 0, fontface = "plain", size = text_size / size_ratio) +
          guides(fill = guide_colorbar(title.vjust = .8)) +
          xlab("RT (min)") +
          ylab("HI (%ACN)")
      }
    }
    
    # Add legend for samples with fake points (outside graph)
    # TODO : Make sure colors defined here correspond exactly with analysis colors later (seems ok for now)
    if (scope == "experiment" && graph_type != "shift_distrib" && pep_nb > 100 * n_analyses) {
      p <- p +
        geom_point(mapping = aes(x = RT + 2 * max_RT, color = Sample)) +
        xlim(NA, max_RT)
    }

    if (n_analyses > 1) {
      if (n_samples > 1) {
        nrows_ncols <- determine_row_col_nb(n_samples, max_ana_per_samp)
        nrows <- nrows_ncols$nrows
        ncols <- nrows_ncols$ncols
        
        # Divide graphs into multiple facets for each analysis
        p <- p +
          facet_wrap(~ graph_label, ncol = ncols, nrow = nrows, scales = "fixed", labeller = labeller(graph_label = strip_sample)) +
          theme(strip.text.x = element_text(size = text_size, face = "italic", margin = margin()),
                strip.text.y = element_text(size = text_size, face = "italic", margin = margin()),
                legend.position = "bottom"
          ) +
          guides(color = guide_legend(order = 1, override.aes = list(size=0.8)))

        fig_width_height <- determine_fig_width_height(nrows, ncols, n_samples, n_analyses, max_ana_per_samp)
        fig_width <- fig_width_height$fig_width
        fig_height <- fig_width_height$fig_height
      } else {
        p <- p +
          facet_wrap(~ Analysis, scales = "fixed", labeller = label_value) +
          theme(strip.text.x = element_text(size = text_size, face = "italic", margin = margin()),
                legend.position = "bottom"
          ) +
          guides(color = guide_legend(order = 1, override.aes = list(size=0.8)))

        fig_width_height <- determine_fig_width_height(NULL, NULL, n_analyses, n_samples, max_ana_per_samp)
        fig_width <- fig_width_height$fig_width
        fig_height <- fig_width_height$fig_height
      }
    } else {
      if (scope == "analysis" && graph_type == "shift_distrib") {  # Do not display sample legend if analysis scope
        p <- p +
          facet_wrap(~ Analysis, scales = "fixed", labeller = label_value) +
          theme(strip.text.x = element_text(size = text_size, face = "italic", margin = margin()),
                legend.position = "none"
          )
      } else {
        p <- p +
          facet_wrap(~ Analysis, scales = "fixed", labeller = label_value) +
          theme(strip.text.x = element_text(size = text_size, face = "italic", margin = margin()),
                legend.position = "bottom"
          ) +
          guides(color = guide_legend(order = 1, override.aes = list(size=0.8)))
      }
    }

    # Change color of facets labels background
    if (scope == "experiment" && graph_type != "shift_distrib" && pep_nb > 100 * n_analyses) {
      p <- p +
        guides(color = guide_legend(title.position = "top", 
                                    title.hjust = 0.5, 
                                    order = 1, 
                                    override.aes = list(size = 0.8)
                                    )
        ) + 
        theme(legend.text = element_text(size = 7))
      ana_colors <- df %>% select(c("Analysis", "Sample", "samp_color")) %>% distinct()
      plot_grob <- ggplotGrob(p)
      plot_grob$grobs[which(grepl("strip-t", plot_grob$layout$name))] <- lapply(plot_grob$grobs[which(grepl("strip-t", plot_grob$layout$name))], change_facet_label_color, ana_colors)
      ggsave(plot_file, plot_grob, device = plot_as, width = fig_width, height = fig_height, units = "cm")
    } else {
      ggsave(plot_file, p, device = plot_as, width = fig_width, height = fig_height, units = "cm")
    }  
  }

  # Depending on the scope, keep only the relevant rows, then draw
  filter_and_draw <- function(df, graph_type, lm_df, scope, plot_files, id) {
    if (scope == "sample") {
      samp_df <- df %>% filter(Sample_ID == id)
      samp_lm_df <- lm_df %>% filter(Sample_ID == id)
      draw_graph(samp_df, graph_type, samp_lm_df, plot_files[[as.character(id)]], scope)
    } else if (scope == "analysis") {
      ana_df <- df %>% filter(Analysis_ID == id)
      ana_lm_df <- lm_df %>% filter(Analysis_ID == id)
      draw_graph(ana_df, graph_type, ana_lm_df, plot_files[[as.character(id)]], scope)
    }
  }

  # Actually do everything here in a kind of vectorized way   
  # Experiment scope
  exp_plot_file <- paste0(plot_file, "_exp.", plot_as)
  draw_graph(df, "HI_vs_RT", lm_df, exp_plot_file, "experiment")
  exp_plot_file_shift <- paste0(plot_file, "_exp_shift.", plot_as)
  draw_graph(df, "shift_distrib", lm_df, exp_plot_file_shift, "experiment")
  exp_plot_file_shift_RT <- paste0(plot_file, "_exp_shift_RT.", plot_as)
  draw_graph(df, "shift_vs_RT", lm_df, exp_plot_file_shift_RT, "experiment")
  
  # Sample scope
  sampIDs <- unlist(distinct(df %>% select(Sample_ID)))
  
  samp_plot_files <- lapply(sampIDs, function(x) paste0(plot_file, "_samp_", x, ".", plot_as))
  names(samp_plot_files) <- sampIDs
  invisible(lapply(sampIDs, function(x) filter_and_draw(df, "HI_vs_RT", lm_df, "sample", samp_plot_files, x)))
  
  samp_plot_files_shift <- lapply(sampIDs, function(x) paste0(plot_file, "_samp_shift_", x, ".", plot_as))
  names(samp_plot_files_shift) <- sampIDs
  invisible(lapply(sampIDs, function(x) filter_and_draw(df, "shift_distrib", lm_df, "sample", samp_plot_files_shift, x)))
  
  samp_plot_files_shift_RT <- lapply(sampIDs, function(x) paste0(plot_file, "_samp_shift_RT_", x, ".", plot_as))
  names(samp_plot_files_shift_RT) <- sampIDs
  invisible(lapply(sampIDs, function(x) filter_and_draw(df, "shift_vs_RT", lm_df, "sample", samp_plot_files_shift_RT, x)))

  
  # Analysis scope
  anaIDs <- unlist(distinct(df %>% select(Analysis_ID)))

  ana_plot_files <- lapply(anaIDs, function(x) paste0(plot_file, "_ana_", x, ".", plot_as))
  names(ana_plot_files) <- anaIDs
  invisible(lapply(anaIDs, function(x) filter_and_draw(df, "HI_vs_RT", lm_df, "analysis", ana_plot_files, x)))

  ana_plot_files_shift <- lapply(anaIDs, function(x) paste0(plot_file, "_ana_shift_", x, ".", plot_as))
  names(ana_plot_files_shift) <- anaIDs
  invisible(lapply(anaIDs, function(x) filter_and_draw(df, "shift_distrib", lm_df, "analysis", ana_plot_files_shift, x)))

  ana_plot_files_shift_RT <- lapply(anaIDs, function(x) paste0(plot_file, "_ana_shift_RT_", x, ".", plot_as))
  names(ana_plot_files_shift_RT) <- anaIDs
  invisible(lapply(anaIDs, function(x) filter_and_draw(df, "shift_vs_RT", lm_df, "analysis", ana_plot_files_shift_RT, x)))
}


########################
### Main computation ###
########################
# Command line arguments and setup
cmd_args <- commandArgs(trailingOnly = TRUE)
work_dir <- cmd_args[1]  # The directory where to find input and output files
info_file <- file.path(work_dir, "hydro_info.txt")
setwd(work_dir)
print(sessionInfo())

if (!file.exists(info_file)) {
  stop(paste("Information file", info_file, "does not exist, make sure that the right path was provided and that it has been written.", sep = " "))
}
params <- read_params(info_file)

hydro_df <- compute_hydro(params$data_file, compute_HI_v1)
lin_mod <- compute_linear_models(hydro_df, params$rt_offset)
hydro_df <- compute_percentile_time_shift(hydro_df, lin_mod)
plot_HI_RT_shifts(hydro_df, lin_mod, params$graph_file, params$plot_as, params$rt_offset)


hydro_df <- hydro_df %>%
  select(-Analysis, -Sample)
lin_mod <- lin_mod %>%
  select(-Analysis, -Sample, -label)

write.table(hydro_df, params$output_hydro, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(lin_mod, params$output_lin_mod, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


####>Revision history<####
# 1.0.0 Created (VL 25/01/2021)
