################################################################################
# GSEA.R       1.0.1                                                           #
# Authors: Victor Laigle (Institut Curie)                                      #
# Contact: myproms@curie.fr                                                    #
# GSEA computation with clusterProfiler package :                              #
# Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for         #
# comparing biological themes among gene clusters.” OMICS: A Journal of        #
# Integrative Biology, 16(5), 284-287. doi: 10.1089/omi.2011.0118              #
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

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

select <- dplyr::select  # potential conflict with other packages
sessionInfo()

cmd_args <- commandArgs(trailingOnly = TRUE)
gsea_dir <- cmd_args[1]
param_file <- file.path(gsea_dir, "gsea_info.txt")
setwd(gsea_dir)


read_params <- function(file) {
  # Take the gsea_info.txt file with parameters to run GSEA and read it.
  # It must be a tsv file with parameter name/value pairs.

  info_lines    <- readLines(file)
  param_names   <- character(length(info_lines))
  param_values  <- character(length(info_lines))

  for (i in 1:length(info_lines)) {
    info  <- strsplit(info_lines[i], split = '\t')[[1]]
    param <- info[1]
    value <- info[-1]

    param_names[i]  <- param
    param_values[i] <- value
  }
  params <- as.list(param_values)
  names(params) <- param_names

  params$is_ratio         <- as.logical(as.integer(params$is_ratio))
  params$pep_weights      <- as.logical(as.integer(params$pep_weights))
  params$pval_cutoff      <- as.double(params$pval_cutoff)
  params$max_relevant_gs  <- as.integer(params$max_relevant_gs)

  return(params)
}


load_data <- function(file, genes_file, weight = FALSE) {
  # Load quantitative data from file
  # Remove duplicate proteins / gene names (typically isoforms)
  # and sort gene_list by decreasing order of quantitative values
  # e.g. ratio of differential expression or absolute amounts
  expression_data <- read_tsv(file, col_types = "ccddd")  # Protein_ID, Gene_name, Value, Peptides_nb
  if (weight) {
    expression_data <- transform_infinite(expression_data)
    expression_data <- compute_weighted_scores(expression_data)
  } else {
    expression_data <- remove_infinite(expression_data)
  }
  expression_data <- remove_duplicates(expression_data)  # Duplicate gene names

  write.table(expression_data, genes_file, sep = "\t", row.names = FALSE, quote = FALSE)

  gene_list <- expression_data$Value
  names(gene_list) <- expression_data$Gene_name
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}


remove_duplicates <- function(data) {
  # Remove duplicate gene names from the matrix
  # The protein kept is the one with the most peptides
  # On equal number of peptides, we keep the most extreme value
  data_dup <- data
  data_no_dup <- data_dup %>%
    arrange(Gene_name, desc(Peptides_nb), desc(abs(Value))) %>%
    distinct(Gene_name, .keep_all = TRUE)

  return(data_no_dup)
}


transform_infinite <- function(data) {
  # Attribute value of the 99th percentile to infinite ratios
  # We consider the absolute log2 values of both up and down ratios
  my_data <- data
  tol <- 1e-8  # Tolerance for numeric comparison
  infinite_ratio <- 1000
  log_infinite = log2(infinite_ratio)

  my_data <- my_data %>%
    mutate(is_inf = (abs(abs(Value) - log_infinite) < tol))  # Not "==" because of floating point issues

  value_99 <- my_data %>%
    filter(!is_inf) %>%
    summarise(percent_99 = quantile(abs(Value), probs = seq(0, 1, by= 0.01))["99%"]) %>%
    pull() %>%
    unname()

  my_data <- my_data %>%
    mutate(Value = case_when(is_inf ~ sign(Value) * value_99, TRUE ~ Value))

  return(my_data)
}


compute_weighted_scores <- function(data) {
  # Weight data according to peptide numbers
  weighted_data <- data

  weighted_data <- weighted_data %>%
    # mutate(Value = Value * Peptides_nb)  # Scoring seems very biased towards prots with high nb of peptides
    mutate(Value = Value * log10(Peptides_nb + 1))  # Other scoring possibility (+1 to keep prots with only 1 pep)
    # mutate(Value = Value * max(10, Peptides_nb))  # Other scoring possibility
    # mutate(Value = Value * -log10(Pvalue)  # Other scoring possibility would be to use the pvalue from DE analysis

  return(weighted_data)
}


remove_infinite <- function(data) {
  # Remove infinite protein ratios before computation of GSEA.
  # Used when protein scores are not weighted
  data_no_inf <- data
  tol <- 1e-8  # Tolerance for numeric comparison
  infinite_ratio = 1000
  log_infinite = log2(infinite_ratio)

  data_no_inf <- data_no_inf %>%
    filter(!(abs(abs(Value) - log_infinite) < tol))  # Not "==" because of floating point issues

  return(data_no_inf)
}


fetch_gene_sets <- function(gmt_file = NULL, gmt_type = NULL) {
  gene_sets <- NULL
  term_gene_map <- NULL
  term_name_map <- NA

  if (is.null(gmt_file)) {
    stop("Cannot get gene sets if no file containing the custom gene sets is provided !")
  } else {
    if (!is.null(gmt_type) && tolower(gmt_type) == "wikipathways") {
      wp2gene <- as_tibble(read.gmt(gmt_file))
      wp2gene <- wp2gene %>% separate(ont, c("name","version","wpid","org"), "%")
      term_gene_map <- wp2gene %>% select(wpid, gene) %>% unique()
      term_name_map <- wp2gene %>% select(wpid, name) %>% unique()
    } else {
      gene_sets <- as_tibble(read.gmt(gmt_file))  # cols: ont, gene (term, gene in latest versions of clusterProfiler)
      term_gene_map <- gene_sets %>% select(ont, gene)  # select(term, gene)
      # term_name_map <- gene_sets %>% select(gs_id, gs_name)
    }
  }
  if (is.null(term_gene_map) || nrow(term_gene_map) == 0) {
    stop("Could not retrieve the gene sets from the specified gene sets databank or file.")
  }
  return(list(t2g = term_gene_map, t2n = term_name_map))
}


plot_gsea <- function(gsea_data, max_gs, out_graph) {
  max_plot <- min(max_gs, nrow(gsea_data))

  for (i in 1:max_plot) {  # GSEA results sorted by p-values automatically (not adjusted p-value)
    geneSetName <- str_to_title(tolower(gsub("_", " ", gsea_data$ID[i], perl = TRUE)))
    gsea_plot <- gseaplot2(gsea_data, geneSetID = i, title = geneSetName, color = "#00BFC4", ES_geom = "line")
    graph_name <- paste0(out_graph, "_", gsea_data$ID[i], ".png")
    ggsave(graph_name, plot = gsea_plot)
  }

  d_plot <- dotplot(gsea_data,
    x = "NES",
    color = "p.adjust",
    showCategory = max_plot,
    size = "setSize",
    split = 2,
    font.size = 9,
    title = "Normalised Enrichment Score by Gene Set"
  )
  d_plot <- d_plot + geom_vline(xintercept = 0, linetype="dashed", color = "darkorange")
  dplot_name <- paste0(out_graph, "_dotplot.png")
  ggsave(dplot_name, plot = d_plot, width = 8, height = 8)
}


########################
### Main Computation ###
########################
if (!file.exists(param_file)) {
  stop(paste("Information file", param_file, "does not exist, make sure that the right path was provided and that it has been written.", sep = " "))
}
params <- read_params(param_file)
# print(params)  # DEBUG

# Load quantification results (ratios etc.)
gene_scores <- load_data(params$data_file, params$genes_used_file, params$pep_weights)

# Load Gene Sets
gene_sets <- fetch_gene_sets(gmt_file = params$gmt_file)
term_gene_map <- gene_sets$t2g
term_name_map <- gene_sets$t2n

# Perform GSEA
gsea_data <- GSEA(
  geneList = gene_scores,
  exponent = 1,
  nPerm = 1000000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = params$pval_cutoff,
  pAdjustMethod = params$pval_correction,
  TERM2GENE = term_gene_map,
  TERM2NAME = term_name_map,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)

if (nrow(gsea_data@result) == 0) {  # No results
  to_write <- c("Not a single Gene Set enriched under the specified p-value cutoff.")
  fileConn <- file(params$output_file)
  writeLines(to_write, fileConn)
  close(fileConn)
} else {  # GSEA results are ok
  tbl_gsea <- as_tibble(as.data.frame(gsea_data))
  write.table(tbl_gsea, params$output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  plot_gsea(gsea_data, params$max_relevant_gs, params$output_graph)
}


### Possible gene sets with clusterProfiler
# MisgDB
# GO Biological Processes
# GO Cellular Components
# GO Molecular Function
# KEGG pathways
# KEGG modules
# WikiPathways
# Cell Marker
# Reactome
# MeSH
# Disease Ontology (DO)
# Network of Cancer Gene (NCG)
# Disease Gene Network (DGN)
# DAVID

####>Revision history<####
# 1.0.1 [MODIF] Add score values to genes_used_file (VL 31/05/21)
# 1.0.0 Created (VL 19/10/20)
