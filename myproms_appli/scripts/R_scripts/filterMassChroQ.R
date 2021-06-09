#!/usr/bin/env Rscript

################################################################################
# filterMassChroQ.R         1.1.0                                              #
# Authors: Victor Laigle (Institut Curie)                                      #
# Contact: myproms@curie.fr                                                    #
# Filter MassChroQ results to remove peptides reextracted with low confidence  #
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

library(methods)
library(ggplot2)
library(scales)  # For hue_pal()
library(dplyr)
library(tidyr)
library(purrr)
library(tictoc)  # To record timing

tic.clearlog()
tic("Whole process")
# User parameters
cmd_args <- commandArgs(trailingOnly = TRUE)
quantif_dir <- cmd_args[1]
setwd(file.path(quantif_dir, "result_quanti.d"))

# File names
input_file <- "peptides_q1_G1.tsv"
output_file <- "peptides_q1_G1_filtered.tsv"
plot_file <- "peptides_rt_sd_hist.png"
filter_info_file <- "filtering_info.txt"
log_file <- "filtering_log.out"

write("1. Processing filter parameters", log_file, append = TRUE)

method <- cmd_args[2]  # One of : no_filter/fixed/max_2FWHM/max_3FWHM/max_4FWHM
limit_threshold <- cmd_args[3]  # Number, considered as NULL if omitted or anything else than a number
min_threshold <- cmd_args[4]  # Number, considered as NULL if omitted or anything else than a number

# Setup for command line arguments
method_choices <- c("no_filter", "fixed", "max_2FWHM", "max_3FWHM", "max_4FWHM")
if (!(method %in% method_choices)) {  # Set method as no_filter if not in the possibilities
  method <- "no_filter"
}

if (is.na(as.numeric(limit_threshold))) {
  limit_threshold <- NULL
} else {
  limit_threshold <- as.numeric(limit_threshold)
}

if (is.na(as.numeric(min_threshold))) {
  min_threshold <- NULL
} else {
  min_threshold <- as.numeric(min_threshold)
}

# Outlier filtering parameters, change here if needed
out_type <- "extrem"  # "extrem" = 3*IQR or "mild" = 1.5*IQR
out_min_sd_increase <- 0.1  # Number between 0.0 and 1.0 (percentage)


find_outliers <- function(full_data, outlier_type = "extrem", min_sd_increase = 0.1) {
  # For each peptide found in multiple samples, find whether there are any outlier
  # (based on retention time) across the different samples.
  # If the peptide has been identified in three samples or less, we consider that there are no outlier.
  # We use the Tukey's fence method to find potential outliers.
  # A peptide is considered as a mild outlier if its retention time is not between Q1 - 1.5 * IQR and Q3 + 1.5 * IQR
  # A peptide is considered as an extreme outlier if its retention time is not between Q1 - 3 * IQR and Q3 + 3 * IQR
  # The effect of each potential outlier on the retention time standard deviation is then assessed
  # and the outlier is removed only if we observe a strong effect on the sd.
  # Finally, a peptide is considered an outlier if :
  #   - it is outside the Tukey's fences on a boxplot (mild or extrem outlier), and
  #   - the change of standard deviation on retention times is greater than the min_sd_increase defined (default 10%)
  #     when comparing the sd with all peptides to the sd with the potential outlier removed.

  if (outlier_type == "mild") {  # Filter on mild outliers
    IQR_weight <- 1.5
  } else {  # Filter only extrem outliers as soon as "mild" was not specifically required
    IQR_weight <- 3
  }

  sd_all_but_one <- function(n, full_grp_data) {  # Computes sd for the subset with peptide n removed
    return(sd(full_grp_data$rt[full_grp_data$N != n], na.rm=TRUE))
  }

  find_outlier_in_grp <- function(grp_data, ...) {
    # Actual function implementing the search for outliers described above, for a given peptide in all samples
    bxplt <- boxplot(grp_data$rt, range = IQR_weight, plot = FALSE)
    sd_all <- sd(grp_data$rt, na.rm = TRUE)

    grp_data <- grp_data %>%
      mutate(potential_outliers = rt %in% bxplt$out)  # Peptides outside the Tukey's fences considered as potential outlier
    
    if (nrow(grp_data) <= 3 || !any(grp_data$potential_outliers)) {  # Don't consider outliers if less than 3 peptides across samples
      grp_data <- grp_data %>%
        mutate(sd_without = NA, sd = sd_all, sd_increase = NA, outlier = FALSE)
    } else {
      grp_data <- grp_data %>%
        mutate(N = 1:n()) %>%
        mutate(sd_without = map(.$N, sd_all_but_one, .)) %>%
        unnest() %>%  # Works with dplyr version of image myproms_1.3.2.img
        # unnest(cols = c(sd_without))  # For newer dplyr versions. Does not have the expected behavior on myproms_1.3.2.img
        mutate(sd = sd_all, sd_increase = sd_all / sd_without - 1, N = NULL)

      grp_data <- grp_data %>%
        mutate(outlier = potential_outliers & (sd_increase > min_sd_increase),
               potential_outliers = NULL)
    }

    return(grp_data)
  }

  data_with_outliers <- full_data %>%
    group_by(peptidez) %>%
    nest() %>%
    mutate(., data = map(data, find_outlier_in_grp)) %>%
    unnest()  # unnest(cols = c(data)) when newer versions of dplyr will be used

  return(data_with_outliers)
}


compute_FWHM <- function(ddata) {  # Input is density data, with x and y columns
  xmax <- ddata$x[ddata$y == max(ddata$y)]
  x1 <- ddata$x[ddata$x < xmax][which.min(abs(ddata$y[ddata$x < xmax] - (max(ddata$y) / 2)))]
  x2 <- ddata$x[ddata$x > xmax][which.min(abs(ddata$y[ddata$x > xmax] - (max(ddata$y) / 2)))]
  
  p1 <- c(x1, ddata$y[ddata$x == x1])
  p2 <- c(x2, ddata$y[ddata$x == x2])
  pmax <- c(xmax, ddata$y[ddata$x == xmax])
  FWHM <- x2 - x1
  
  return(list(FWHM = FWHM, start = p1, end = p2, pmax = pmax))
}


compute_threshold <- function(how = "max_2FWHM", FWHM = NULL, max = NULL, limit_thresh = NULL, min_thresh = 10) {
  threshold <- -1
  if (how == "fixed") {
    if (is.null(limit_thresh)) {
      stop("Provide a threshold if you want to set it yourself ! The limit_thresh parameter has to be set.")
    } else {
      threshold <- limit_thresh
    }
  } else if (how == "max_2FWHM") {
    threshold <- max + 2 * FWHM
  } else if (how == "max_3FWHM") {
    threshold <- max + 3 * FWHM
  } else if (how == "max_4FWHM") {
    threshold <- max + 4 * FWHM
  } else if (how == "no_filter") {
    threshold <- -1
  }
  
  if (!is.null(limit_thresh) && threshold > limit_thresh) {
    threshold <- limit_thresh
  } else if (how != "fixed" && how != "no_filter" && !is.null(min_thresh) && threshold < min_thresh) {
    threshold <- min_thresh
  }
  
  return(threshold)
}


get_threshold <- function(rt_data, graph_file, method, limit_threshold, min_threshold, total_nb) {
  ### Density computation
  # Discard peptides where dispersion of retention time is greater than the 95th percentile 
  # for density and threshold calculations. The histograms I've seen that far have a very 
  # long right tail due to a few peptides with very high standard deviation on RT, 
  # which needs to be cut to get proper computation of density. 
  # If the values go too far to the right, the density computation approximates the histogram very poorly,
  # probably because the density function relies on a normal approximation,
  # which is completely lost with a very skewed distribution.
  # 95% seems like a good tradeoff from what I've seen but it is arbitrary. 
  # It still allows to compute the density on the important part (the peak) and the FWHM.
  percentile_95 <- quantile(rt_data$rt_sd, probs = seq(0, 1, by= 0.01))["95%"]
  percentile_1 <- quantile(rt_data$rt_sd, probs = seq(0, 1, by= 0.01))["1%"]
  rt_data_1_95 <- rt_data[rt_data$rt_sd > percentile_1 & rt_data$rt_sd <= percentile_95, ]
  
  hist_bw <- max(rt_data_1_95$rt_sd) / 100
  rt_sd_density <- density(rt_data_1_95$rt_sd, bw = "SJ", adjust = 1.5)
  count_scale <- hist_bw * nrow(rt_data_1_95)  # To scale from density to counts
  rt_sd_density.data <- data.frame(x = rt_sd_density$x, y = rt_sd_density$y * count_scale)
  
  ## Computation of FWHM from density data
  FWHM_data <- compute_FWHM(rt_sd_density.data)
  # Need creation of data.frame here to get ggplot and ggsave to work (cf. data parameter)
  FWHM.data <- data.frame(FWHM = FWHM_data$FWHM,
                          start_x = FWHM_data$start[1],
                          start_y = FWHM_data$start[2],
                          max_x = FWHM_data$pmax[1],
                          max_y = FWHM_data$pmax[2],
                          end_x = FWHM_data$end[1],
                          end_y = FWHM_data$end[2]
  )

  threshold <- compute_threshold(how = method, 
                                 FWHM = FWHM.data$FWHM,
                                 max = FWHM.data$max_x,
                                 limit_thresh = limit_threshold,
                                 min_thresh = min_threshold
  )

  FWHM.data$threshold <- threshold

  ## Plot the histogram of peptides RT dispersion values
  hist_rt_sd <- NULL
  if (threshold < 0) {  # No threshold and no filter on peptides RT dispersion, just display histogram
    # Plot histogram of peptides retention times dispersion
    hist_rt_sd <- ggplot(rt_data, aes(x = rt_sd)) +
      geom_histogram(aes(y = ..count.., color = hue_pal()(1), fill = hue_pal()(1)), binwidth = hist_bw, alpha = 0.7, na.rm=TRUE) +
      xlim(c(0, 200)) +
      labs(title = "Dispersion of peptides retention times", x = "RT Standard deviation", y = "Counts") + 
      theme(legend.position="none")
  } else {
    kept_nb <- nrow(rt_data[rt_data$rt_sd <= threshold, ])
    filtered_nb <- total_nb - kept_nb
    filtered_pct <- round(filtered_nb / total_nb * 100, 1)

    graph_xlim <- ifelse(threshold <= 100, ifelse(threshold <= 40, 50, 100), round(threshold + 20, -1))
    
    # Plot histogram of peptides retention times dispersion + other info
    # TODO : aggregate annotate "text" into a geom_text element (reduce code redundancy)
    hist_rt_sd <- ggplot(rt_data, aes(x = rt_sd)) +
      geom_histogram(aes(y = ..count.., color = hue_pal()(1), fill = hue_pal()(1)), binwidth = hist_bw, alpha = 0.7, na.rm=TRUE) +
      geom_line(aes(x = x, y = y), rt_sd_density.data, alpha = 0.8, color="darkred", na.rm=TRUE) +
      geom_segment(aes(x = start_x,
                       y = mean(start_y, end_y),
                       xend = end_x,
                       yend = mean(start_y, end_y)),
                   FWHM.data,
                   color = "darkred",
                   size = 0.2,
                   na.rm=TRUE
      ) +
      geom_vline(aes(xintercept = threshold), FWHM.data, color = "darkblue", na.rm=TRUE) +
      xlim(c(0, graph_xlim)) + 
      labs(title = "Threshold on the dispersion of peptides retention times", x = "RT Standard deviation", y = "Counts") + 
      theme(legend.position="none") +
      annotate("text",
               x = FWHM.data$threshold + 1,
               y = FWHM.data$max_y,
               label = paste0("Threshold : ", format(round(FWHM.data$threshold, 2), nsmall = 2), " sec\n",
                              "Peptides accepted : ", kept_nb,
                              "\nOut of ", total_nb, " (distinct peptides across all analyses)\n",
                              filtered_nb, " (", filtered_pct, "%) peptides filtered, besides outliers."),
               size = 5,
               fontface = "italic",
               hjust = 0,  # Horizontal justification = left
               vjust = 1  # Vertical justfication = top
      ) +
      annotate("text",
               x = FWHM.data$end_x + 1,
               y = FWHM.data$end_y,
               label = paste0("FWHM\n", format(round(FWHM.data$FWHM, 2), nsmall = 2), " sec"),
               size = 4,
               fontface = "italic",
               hjust = 0  # Horizontal justification = left
      )
  }
  suppressWarnings(ggsave(graph_file, plot = hist_rt_sd, width = 10, height = 10))
  
  return(list(thresh = threshold, FWHM = FWHM.data$FWHM))
}


### Start main computation ###
write("2. Reading MassChroQ data", log_file, append = TRUE)
mcq_out <- read.table(file = input_file, sep = '\t', header = TRUE)
mcq_filtered <- NULL
filter_info_param <- NULL
filter_info_value <- NULL

if (nrow(unique(mcq_out[c("quantification", "group", "msrun", "msrunfile")])) > 1) {  # Check if at least two samples were aligned, otherwise it is impossible to compute std deviation
  initial_pep_nb <- nrow(mcq_out)
  mcq_out$peptidez <- paste(mcq_out$peptide, mcq_out$z, sep="-")
  keep_cols <- colnames(mcq_out)

  write("3/6 Checking and removing outliers (if needed)", log_file, append = TRUE)
  nb_non_outliers <- NULL
  if (method == "no_filter") {  # Remove outliers only if filtering
    nb_non_outliers <- initial_pep_nb
  } else {
    outliers_data <- find_outliers(mcq_out, outlier_type = out_type, min_sd_increase = out_min_sd_increase)
    non_outliers <- outliers_data[!outliers_data$outlier, ]
    nb_non_outliers <- nrow(non_outliers)
  
    mcq_out <- inner_join(mcq_out, non_outliers, by = keep_cols)[, keep_cols]  # Keep only the non-outliers
  }
  outlier_pep_nb <- initial_pep_nb - nb_non_outliers
  keep_cols <- keep_cols[keep_cols != "peptidez"]
  mcq_filtered <- mcq_out[, keep_cols]

  total_rt <- aggregate(mcq_out$rt, list(peptidez = mcq_out$peptidez), FUN=sd)
  pep_nb <- nrow(total_rt)
  total_rt <- na.omit(total_rt)
  colnames(total_rt) <- c("peptidez", "rt_sd")

  write("4/6 Computing FWHM and threshold on RT standard deviation", log_file, append = TRUE)
  filter_data <- get_threshold(total_rt, plot_file, method = method, limit_threshold = limit_threshold, min_threshold = min_threshold, total_nb = pep_nb)
  threshold <- filter_data$thresh
  FWHM <- filter_data$FWHM

  # Discard peptides having a retention time standard deviation greater than the computed threshold (e.g. 15s)
  write("5/6 Removing peptides which do not pass the filter", log_file, append = TRUE)
  nb_unique_after_filter <- NA
  if (threshold >= 0) {
    peptidez_with_good_rt <- total_rt$peptidez[total_rt$rt_sd <= threshold]
    nb_unique_after_filter <- length(peptidez_with_good_rt)
    # Keep peptides with SD < threshold but also those with SD = NA (not in total_rt = only one sample had this peptide)
    mcq_filtered <- mcq_filtered[(mcq_out$peptidez %in% peptidez_with_good_rt) | !(mcq_out$peptidez %in% total_rt$peptidez), ]
  } else {
    nb_unique_after_filter <- pep_nb
  }
  final_pep_nb <- nrow(mcq_filtered)
  pep_filtered_thresh <- nb_non_outliers - final_pep_nb
  total_pep_filtered <- initial_pep_nb - final_pep_nb

  # Writing filtering info
  write("6/6 Writing filtered peptides table and filtering info", log_file, append = TRUE)
  filter_info_param <- NULL
  filter_info_value <- NULL
  if (method == "no_filter") {  # Reduced info since filter was not applied
    filter_info_param <- c(
      "Initial number of peptide ions",
      "Unique peptide ions across all samples",
      "Filtering method",
      "Full Width at Half Maximum of RT SD distribution (s)"
    )
    filter_info_value <- c(
      initial_pep_nb,
      pep_nb,
      method,
      format(round(FWHM, 2), nsmall = 2)
    )
  } else {
    if (is.null(limit_threshold)) {  # The filtering info is shifted if we try to write a NULL parameter
      limit_threshold <- "None"
    }
    if (is.null(min_threshold)) {  # Filtering info shifted when trying to write NULL parameter
      min_threshold <- "None"
    }
    filter_info_param <- c(
      "Initial number of peptide ions",
      "Unique peptides ions across all samples",
      "Filtering method",
      "Minimum threshold on RT standard deviation (s)",
      "Maximum threshold on RT standard deviation (s)",
      "Type of outliers (for Tukey fences)",
      "Minimum SD increase to confirm outlier",
      "Number of outlier peptide ions",
      "Full Width at Half Maximum of RT SD distribution (s)",
      "Computed threshold (s)",
      "Peptide ions filtered at computed threshold",
      "Total number of peptide ions filtered (outliers + threshold)",
      "Valid peptide ions after filtering",
      "Valid unique peptide ions (across all samples) after filtering"
    )
    filter_info_value <- c(
      initial_pep_nb,
      pep_nb,
      method,
      min_threshold,
      limit_threshold,
      out_type,
      paste0(round(out_min_sd_increase * 100, 1), "%"),
      outlier_pep_nb,
      format(round(FWHM, 2), nsmall = 2),
      format(round(threshold, 2), nsmall = 2),
      pep_filtered_thresh,
      total_pep_filtered,
      final_pep_nb,
      nb_unique_after_filter
    )
    if (method == "fixed") {
      filter_info_param[[5]] <- "Fixed threshold on RT standard deviation (s)"
    }
  }
} else {  # Only one sample : just rewrite MassChroQ output as is + minimal filtering info
  mcq_filtered <- mcq_out
  initial_pep_nb <- nrow(mcq_out)

  write("3/3 Writing peptides table and info", log_file, append = TRUE)
  filter_info_param <- c(
    "Total number of peptide ions in analysis",
    "Filtering method"
  )
  filter_info_value <- c(
    initial_pep_nb,
    method
  )
}

write.table(mcq_filtered, file = output_file, quote = FALSE, sep = "\t", na = "", row.names = FALSE, col.names = TRUE)

filtering_info <- paste(filter_info_param, filter_info_value, sep = "\t")
writeLines(filtering_info, filter_info_file)

toc(log = TRUE, quiet = TRUE)
tictoc_sec <- round(tic.log(format = FALSE)[[1]]$toc - tic.log(format = FALSE)[[1]]$tic, 2)
tictoc_min <- round(tictoc_sec / 60, 1)
write(paste0("Ended. The whole process took ", tictoc_sec, "s (", tictoc_min, "min)."), log_file, append = TRUE)

####>Revision history<####
# 1.1.0 [ENHANCEMENT] Add minimum threshold to filter + improve computation of density (VL 15/04/21)
# 1.0.9 [BUGFIX] Minor fix in parameter list for no filtering (PP 12/02/21)
# 1.0.8 [ENHANCEMENT] Add logs at each step + [BUGFIX] Avoid wrong shift in filtering info (VL 01/02/21)
# 1.0.7 [MINOR] Change filtering_info strings for display (VL 05/10/20)
# 1.0.6 [BUGFIX] Handle case when only one analysis selected for MassChroQ reextraction with RT filter (VL 28/09/20)
# 1.0.5 [BUGFIX] Change set of peptides on which density and FWHM are computed to avoid wrong density computation (VL 25/09/20)
# 1.0.4 [BUGFIX] Add condition for outliers removal, only when the results are filtered (VL 11/09/20)
# 1.0.3 [MINOR] Change possible filters to force a reference to the peak maximum when filter is computed (VL 27/08/20)
# 1.0.2 Add writing of filtering info to file (VL 17/08/20)
# 1.0.1 Add filtering of outliers prior to sd computation (VL 31/07/20)
# 1.0.0 Created (VL 15/07/2020)
