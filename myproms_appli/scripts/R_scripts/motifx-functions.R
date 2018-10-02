################################################################################
# motifx-functions.R       1.0.0                                               #
# Authors: Stephane Liva    												   #
# Contact: myproms@curie.fr                                                    #
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
# motif-x implementation in R
# Created on : 09-11-14 
# Author     : omarwagih
# Wagih O, Sugiyama N, Ishihama Y, Beltrao P. (2015) Uncovering phosphorylation-based specificities through functional interaction networks (2015). Mol. Cell. Proteomics

##############################
#
#	MOTIF X FUNCTIONS
#
##############################

# All amino acids
.GetAA <- function(){
  AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  AA
}

# Constructing amino acids
.BuildPWM <- function(seqs, relative.freq=T){
  
  # Ensure same length characters 
  num.pos = seq.len = nchar(seqs[1])
  if(F){
    # Slows things down for large number of sequences
    seq.len = sapply(seqs, nchar)
    num.pos = seq.len[1]
    if(! all(seq.len == num.pos)) stop('Unequal length of sequences')
  }
  
  # List of valid amino acids, sorted
  namespace = .GetAA()
  
  # Make matrix of letters
  pwm.matrix = sapply(1:num.pos, function(col){
    i = col
    # Get frequencies 
    t = table(substr(seqs, col,col))
    # Match to aa
    ind = match(namespace, names(t))
    # Create column
    col = t[ind]
    col[is.na(col)] = 0
    names(col) = namespace
    
    
    # Do relative frequencies
    if(relative.freq)
      col = col / sum(col)
    
    col
  })
  colnames(pwm.matrix) = 1:num.pos
  
  return(pwm.matrix)
}

.MotifRegex <- function(motifs, cent.regex='[ST]', n=15){
  if(class(motifs) == 'data.frame') motifs = list(motifs)
  
  # Central residue position
  cent.ind = ceiling(n/2)
  # Must have at least one row
  motifs = motifs[sapply(motifs, nrow) > 0]
  regex = sapply(motifs, function(mt){
    aa = mt[[1]]
    col = mt[[2]]
    z = rep('.', n)
    z[col] = aa
    z[cent.ind] = cent.regex
    paste0(z, collapse='')
  })
  regex
}

.FindMotif <- function(fg.seqs, bg.seqs, min.seqs=20, pval.cutoff=0.05, cent.regex='[ST]', verbose=T, perl.impl=F){
  
	AA = .GetAA()
  
  .fg.seqs = fg.seqs
  .bg.seqs = bg.seqs
  
  fg.size = length(fg.seqs)
  bg.size = length(bg.seqs)
  
  kmer.length = nchar(fg.seqs[1])
  cent.index = ceiling(kmer.length/2)
  
  if(verbose) writeLines('Step 1: recursive motif building')
  motif = matrix(0, 0, 2)
  motif.parts = matrix(0, 0, 2)
  pvals = c()
  iter = 1
  while(TRUE){
    # Break if we dont have any sequences remaining
    if(length(bg.seqs) == 0 | length(fg.seqs) == 0) break
    
    # Construct PFM and PWM
    ppwm = .BuildPWM(fg.seqs, relative.freq=F)
    npwm = .BuildPWM(bg.seqs)
    
    # Compute binomial matrix
    bin = 1 - pbinom(ppwm-1, length(fg.seqs), npwm, lower.tail=TRUE)
    
    rownames(bin) = rownames(ppwm)
    colnames(bin) = colnames(ppwm)
    
    # Lowest p-value computed by pbinom is 10e-16, set any zero values to that
    if(perl.impl){
      if(verbose) writeLines('Using perl implementation.')
      bin[bin == 0] = 1e-16 # this is what weblogo uses in perl, different for R
    }else{
      bin[bin == 0] = .Machine$double.xmin # Usually ~= 2.2e-308
    }
    
    
    # Don't enrich on central residue or previous motifs
    bin[,cent.index] = 1
    bin[motif] = 1
    bin[motif.parts] = 1
    
    # bin[, c(1:2, 14:15)] = 1 # exclude outermost
    
    # Give anything with a frequency less than the threshold a p-value of 1
    bin[ppwm < min.seqs] = 1
    
    # Find the row/column of the min binomal p-value
    min.bin = min(bin)
    mbin = which(bin == min.bin & bin < pval.cutoff, arr.ind=T)
    
    # No more significant, breaking
    if(nrow(mbin) == 0) break
    
    # Case where there are more than one matches (likely pvalue=0)
    # Find match with greatest value in PFM
    if(nrow(mbin) > 1){
      r = which.max(ppwm[mbin])
      rn = rownames(mbin)
      mbin = mbin[r,]
      mbin = t(as.matrix(mbin))
      rownames(mbin) = rn[r]
    }
    
    aa = rownames(mbin)[1]
    row = mbin[1,1]
    col = mbin[1,2]
    
    aa. = aa
    
    # Extract sequences 
    ind.pos = substr(fg.seqs, col, col) %in% aa.
    ind.neg = substr(bg.seqs, col, col) %in% aa.
    
    
    if(verbose)
      writeLines(sprintf('\tIteration %s, aa=%s row=%s col=%s, pval=%s, fg=%s, bg=%s', 
                         iter, aa, row, col, signif(min.bin, 3), length(fg.seqs), length(bg.seqs)))
    
    fg.seqs = fg.seqs[ind.pos]
    bg.seqs = bg.seqs[ind.neg]
    motif = rbind(motif, c(row, col))
    
    
    pvals = append(pvals, min.bin)
    iter = iter + 1
  }
  
  # Motif data: data frame with amino acid and positions for the motif
  motif.data = as.data.frame(motif) 
  names(motif.data) = c('aa', 'col')
  
  AA. = AA
  motif.data$aa = AA.[motif.data$aa]
  
  # If we don't have any data, return NULL so that the parent function can deal with it
  if(nrow(motif.data) == 0) return(NULL)
  
  # Find the regex of the motif given the row/col and aa values
  motif.regex = .MotifRegex(motif.data, cent.regex, kmer.length)
  
  # Compute the score of the motif
  motif.score = ifelse( is.null(pvals), NA, sum( -log10(pvals) ) )
  
  # Compute fg/bg matches
  fg.matches = length(grep(motif.regex, .fg.seqs))
  bg.matches = length(grep(motif.regex, .bg.seqs))
  
  return(list(pos=fg.seqs, motif.data=motif.data, motif.score=motif.score, motif.regex=motif.regex,
              fg.matches=fg.matches, fg.size=fg.size, 
              bg.matches=bg.matches, bg.size=bg.size))
}

###############################################
#
# GGSEQLOGO FUNCTIONS
#
###############################################
### Namespaces
##.AA_NAMESPACE = function() c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
##.DNA_NAMESPACE = function() c('A', 'T', 'G', 'C')
##.RNA_NAMESPACE = function() c('A', 'U', 'G', 'C')
##
### Generate letter matrix from vector of sequences
### 
### @param input vector of sequences
##letterMatrix <- function(input){
##  # Ensure kmers are the same length characters 
##  seq.len = sapply(input, nchar)
##  num_pos = seq.len[1]
##  if(! all(seq.len == num_pos)) stop('Sequences in alignment must have identical lengths')
##  
##  # Construct matrix of letters
##  split = unlist( sapply(input, function(seq){strsplit(seq, '')}) )
##  
##  t( matrix(split, seq.len, length(split)/num_pos) )
##}
##
### Guess sequence type based on letter matrix
### 
### @param sp letters
##guessSeqType <- function(sp){
##  # Ensure we have something
##  if(length( intersect(sp, c(.AA_NAMESPACE(), .DNA_NAMESPACE(),.RNA_NAMESPACE())) ) == 0)
##    stop('Could not get guess seq_type. Please explicitly define sequence type or use "other" with custom namespaces.')
##  
##  dat = setdiff(intersect(sp, .AA_NAMESPACE()), c(.DNA_NAMESPACE(),.RNA_NAMESPACE()))
##  if(length(dat) > 0){
##    return('AA')
##  }else if('U' %in% sp){
##    return('RNA')
##  }
##  return('DNA')
##}
##
##
### Find namespace
### 
### @param letter_mat Matrix of latters
### @param seq_type Sequence type
### @param namespace Alphabet
##findNamespace <- function(letter_mat, seq_type, namespace){
##  
##  # Get all letters in our alignment
##  sp = as.character(letter_mat)
##  
##  # Other namespace
##  if(seq_type == "other"){
##    if(is.null(namespace)) 
##      stop('seq_type of "other" must have a defined namespace')
##    
##    namespace = as.character(namespace)
##    # Get unique
##    namespace = unique( unlist(strsplit(namespace, '')) )
##    
##    # Validate
##    non_alphanumeric = grepl("[^[:alnum:] ]", namespace)
##    if( any( non_alphanumeric ) )
##      stop('All letters in the namespace must be alphanumeric')
##    
##    # Ensure there is something in each column
##    # apply(letter_mat, 2, function(column_letters){
##    #   int = intersect(namespace, column_letters)
##    #   if(length(int) == 0)
##    #     stop('The alignment has no letters in namespace match aligned sequences in at least one column')
##    # })
##    
##  }else{
##    if(!is.null(namespace)) 
##      stop('For custom namespaces please set seq_type to "other"')
##    
##    # Guess sequence type
##    if(seq_type == "auto")
##      seq_type = guessSeqType(sp)
##    
##    # Get predefined namespace
##    namespace = get( sprintf('.%s_NAMESPACE', toupper(seq_type)) )()
##  }
##  
##  return(list(seq_type = toupper(seq_type), 
##              namespace = namespace))
##}
##
### Calcualte bits
###
### @param pwm Position weight matrix
### @param N Number of letters in namespace
### @param Nseqs Number of sequences in PWM
##computeBits <- function(pwm, N=4, Nseqs=NULL){
##  H_i = - apply(pwm, 2, function(col) sum(col * log2(col), na.rm=T))
##  e_n = 0
##  if(!is.null(Nseqs)) e_n = (1/logb(2)) * (N-1)/(2*Nseqs) 
## 
##  R_i = log2(N) - (H_i  + e_n)
##  # Set any negatives to 0
##  R_i = pmax(R_i, 0)
##  return(R_i)
##}
##
### Construct relative frequency matrix
### @param seqs aligned sequences as vector
### @param seq_type sequence type
### @param namespace letters used for matrix construction
### @param keep_letter_mat Keep letter matrix for some height methods
##makePFM <- function(seqs, seq_type='auto', namespace=NULL, keep_letter_mat=F){
##  
##  letter_mat = NA
##  if(is.matrix(seqs)){
##    # Process matrix
##    if(is.null(rownames(seqs))) stop('Matrix must have letters for row names')
##    
##    num_pos = ncol(seqs)
##    
##    # Get namespace
##    ns = findNamespace(rownames(seqs), seq_type, namespace)
##    namespace = ns$namespace
##    seq_type = ns$seq_type
##    
##    nseqs = NULL
##    
##    bg_prob = NA
##    pfm_mat = seqs
##    pfm_mat = apply(pfm_mat, 2, function(x) x / sum(x, na.rm=T))
##    
##    missing_rows = setdiff(namespace, rownames(pfm_mat))
##    
##    if(length(missing_rows) > 0){
##      miss = matrix(rep(0, length(missing_rows) * ncol(pfm_mat)), nrow=length(missing_rows), dimnames = list(missing_rows))
##      pfm_mat = rbind(pfm_mat, miss)
##    }
##    
##    pfm_mat = pfm_mat[namespace,]
##
##  }else{
##    # Process sequences
##    
##    # Number of positions in alignment
##    num_pos = nchar(seqs[1])
##    # Number of sequences
##    nseqs = length(seqs)
##    # Letter matrix
##    letter_mat = letterMatrix(seqs)
##    
##    
##    # Get namespace
##    ns = findNamespace(letter_mat, seq_type, namespace=namespace)
##    namespace = ns$namespace
##    seq_type = ns$seq_type
##    
##    # Construct PWM
##    pfm_mat = apply(letter_mat, 2, function(pos.data){
##      # Get frequencies 
##      t = table(pos.data)
##      # Match to aa
##      ind = match(namespace, names(t))
##      # Create column
##      col = t[ind]
##      col[is.na(col)] = 0
##      names(col) = namespace
##      # Do relative frequencies
##      col = col / sum(col)
##      col
##    })
##    
##    mat = matrix((letter_mat %in% namespace), nrow=nrow(letter_mat))
##    attr(pfm_mat, 'nongapped') = apply(mat, 2, sum)/nseqs
##  }
##  
##  # Number of letters in ns
##  N = length(namespace)
##  
##  # Assign seq type and namespace as attributes
##  attr(pfm_mat, 'seq_type') = seq_type
##  attr(pfm_mat, 'namespace') = namespace
##
##  # Non-gapped columns
##  if(seq_type == 'aa') namespace = c(namespace, 'X', 'B', 'Z')
##
##  # Information content
##  attr(pfm_mat, 'bits') = computeBits(pfm_mat, N, nseqs)
##  
##  # Assign AA names to rows/pos col
##  rownames(pfm_mat) = namespace
##  colnames(pfm_mat) = 1:num_pos
##  
##  if(keep_letter_mat) return(list(letter_mat = letter_mat, pfm=pfm_mat))
##
##  return(pfm_mat)
##}
##
##
##
########################
### Matrix to heights
########################
##
### General function to convert matrix of heights to polygon data frame 
### @param mat matrix of heghts
### @param seq_type sequence type
### @decreasing Sets order of letters, high to low or low to high
##matrix_to_heights <- function(mat, seq_type, decreasing=T){
##  
##  mat[is.infinite(mat)] = 0 
##  
##  if(any(duplicated(rownames(mat)))) stop('Matrix input must have unique row names')
##  
##  dat = lapply(1:ncol(mat), function(i){
##    vals = mat[,i]
##    
##    pos = sort( vals[vals >= 0], decreasing = decreasing)
##    neg = sort(vals[vals < 0], decreasing = !decreasing)
##    #vals = sort(vals, decreasing = T)
##    cs_pos = cumsum( pos )
##    cs_neg = cumsum( neg )
##    
##    df_pos = df_neg = NULL
##    
##    if(length(pos) > 0)
##      df_pos = data.frame(letter=names(pos), position=i,  y0=c(0, cs_pos[-length(cs_pos)]), 
##                          y1=cs_pos, stringsAsFactors = F)
##    
##    if(length(neg) > 0)
##      df_neg = data.frame(letter=names(neg), position=i, y0=cs_neg, y1=c(0, cs_neg[-length(cs_neg)]), 
##                          stringsAsFactors = F)
##    
##    rbind(df_pos, df_neg)
##  })
##
##  dat = do.call(rbind, dat)
##
##  # Adjust y spacing 
##  space_factor = 0.004
##  y_pad = max(dat$y1) * space_factor
##  dat$y0 = dat$y0 + y_pad
##  dat = subset(dat, dat$y1 > dat$y0)
##  
##  # Dummy points to make sure full plot is drawn
##  # Make sure position 1 and n have a dummy empty letter missing
##  dummy = data.frame(letter=dat$letter[1], position=NA, y0=0, y1=0)
##  
##  # Missing first position
##  if(dat$position[1] != 1){
##    dummy$position = 1
##    dat = rbind( dummy, dat )
##  }
##  
##  # Missing last position
##  if(dat$position[nrow(dat)] != ncol(mat)){
##    dummy$position = ncol(mat)
##    dat = rbind( dat, dummy )
##  }
##
##  rownames(dat) = NULL
##  
##  attr(dat, 'seq_type') = seq_type
##  
##  dat
##}
##
##
##
### Shannon entropy method
##bits_method <- function(seqs, decreasing, ...){
##  # Make PFM
##  pfm = makePFM(seqs, ...)
##
##  # Get ic
##  ic = attr(pfm, 'bits')
##  if(all(ic == 0)) ic = ic + 2
##  heights = t(t(pfm) * ic)
##
##  seq_type = attr(pfm, 'seq_type')
##  matrix_to_heights(heights, seq_type, decreasing)
##} 
##
### Probability method
##probability_method <- function(seqs, decreasing, ...){
##  # Make PFM
##  pfm = makePFM(seqs, ...)
##  seq_type = attr(pfm, 'seq_type')
##  matrix_to_heights(pfm, seq_type, decreasing)
##}
##
### Change range of values
##newRange <- function(old_vals, new_min=0, new_max=1){
##  old_min = min(old_vals)
##  old_max = max(old_vals)
##  
##  new_vals = (((old_vals - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min
##  new_vals
##}
##
##
###' List fonts available in ggseqlogo
###' 
###' @param v If true, font names are printed to stderr. Otherwise, font names are returned as a character vector
###' @export
##list_fonts <- function(v=T){
##  
##  fonts = c('helvatica','helvetica_regular','helvetica_bold', 'helvetica_light',
##            'roboto_medium','roboto_bold', 'roboto_regular',
##            'akrobat_bold', 'akrobat_regular', 
##            'roboto_slab_bold', 'roboto_slab_regular', 'roboto_slab_light', 
##            'xkcd_regular')
##  if(!v) return(fonts)
##  message('Available ggseqlogo fonts:')
##  for(f in fonts) message('\t', f)
##}
##
##
### Read font from file if not in global envir.
##get_font <- function(font){
##  
##  GGSEQLOGO_FONT_BASE = getOption('GGSEQLOGO_FONT_BASE')
##  if(is.null(GGSEQLOGO_FONT_BASE)){
##    GGSEQLOGO_FONT_BASE=system.file("fonts", "", package = "ggseqlogo")
##    options(GGSEQLOGO_FONT_BASE=GGSEQLOGO_FONT_BASE)
##  }
##  
##  #all_fonts = c('sf_bold', 'sf_regular', 'ms_bold', 'ms_regular', 'xkcd_regular')
##  font = match.arg(tolower(font), list_fonts(F))
##  font_filename = paste0(font, '.font')
##  font_obj_name = sprintf('.ggseqlogo_font_%s', font)
##  
##  font_obj = getOption(font_obj_name)
##  if(is.null(font_obj)){
##    # Not loaded into global env yet - load it into options
##    font_path = file.path(GGSEQLOGO_FONT_BASE, font_filename)
##    font_obj_list = list( tmp=readRDS(font_path) )
##    names(font_obj_list) = font_obj_name
##    options(font_obj_list)
##    font_obj = font_obj_list[[1]]
##  }
##  
##  # Return font data
##  font_obj
##}
##
##
### Generate height data for logo
##logo_data <- function( seqs, method='bits', stack_width=0.95, 
##                       rev_stack_order=F, font, seq_group=1, 
##                       seq_type = 'auto', namespace=NULL ){
##
##  # Get font 
##  font_df = get_font(font)
##  
##  # TODO
##  # hh = twosamplelogo_method(seqs, seqs_bg, pval_thresh=0.05, seq_type = seq_type, namespace = namespace)
##  
##  # Generate heights based on method
##  if(method == 'bits'){
##    hh = bits_method(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace)
##  }else if(method == 'probability'){
##    hh = probability_method(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace)
##  }else if(method == 'custom'){
##    if(seq_type == 'auto') seq_type = guessSeqType(rownames(seqs))
##    hh = matrix_to_heights(seqs, seq_type, decreasing = rev_stack_order)
##  }else{
##    stop('Invalid method!')
##  }
##  
##  # Merge font df and heights
##  ff = merge(font_df, hh, by = 'letter')
##  # Scale x and ys to new positions
##  x_pad = stack_width/2
##  ff$x = newRange(ff$x, ff$position - x_pad, ff$position + x_pad)
##  ff$y = newRange(ff$y, ff$y0, ff$y1)
##  
##  # Rename columns
##  ff = as.data.frame(ff)[,c('x', 'y', 'letter', 'position', 'order')]
##  ff$seq_group = seq_group
##  
##  # Set sequence type as attribute, to be used downstream
##  attr(ff, 'seq_type') = attr(hh, 'seq_type')
##  
##  # Return data table
##  ff
##}
##
###' ggseqlogo custom theme
###' 
###' @param base_size font base size
###' @param base_family font base family
###' 
###' @export
##theme_logo <- function(base_size=12, base_family=''){
##  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace% 
##    theme(panel.grid = element_blank(), legend.position = 'bottom', 
##          axis.text.x=element_text(colour="black"),
##          axis.text.y=element_text(colour="black"))
##}
##
###' Plots sequence logo as a layer on ggplot 
###' 
###' @param data Character vector of sequences or named list of sequences. All sequences must have same width.
###' @param method Height method, can be one of "bits" or "probability" (default: "bits")
###' @param seq_type Sequence type, can be one of "auto", "aa", "dna", "rna" or "other" 
###' (default: "auto", sequence type is automatically guessed)
###' @param namespace Character vector of single letters to be used for custom namespaces
###' @param font Name of font. See \code{list_fonts} for available fonts.
###' @param stack_width Width of letter stack between 0 and 1 (default: 0.95)
###' @param rev_stack_order If \code{TRUE}, order of letter stack is reversed (default: FALSE)
###' @param col_scheme Color scheme applied to the sequence logo. See \code{list_col_schemes} for available fonts.
###' (default: "auto", color scheme is automatically picked based on \code{seq_type}). 
###' One can also pass custom color scheme objects created with the \code{make_col_scheme} function
###' @param low_col,high_col Colors for low and high ends of the gradient if a quantitative color scheme is used (default: "black" and "yellow").
###' @param na_col Color for letters missing in color scheme (default: "grey20")
###' @param plot If \code{FALSE}, plotting data is returned 
###' @param ... Additional arguments passed to layer params
###' 
###' @export
###' @import ggplot2
###' 
###' @examples
###' # Load sample data
###' data(ggseqlogo_sample)
###' 
###' # Produce single sequence logo using geom_logo
###' p1 = ggplot() + geom_logo(seqs_dna[[1]]) + theme_logo()
###' 
##geom_logo <- function(data = NULL, method='bits', seq_type='auto', namespace=NULL,
##                       font='roboto_medium',stack_width=0.95, rev_stack_order=F, col_scheme = 'auto',
##                      low_col='black', high_col='yellow', na_col='grey20',
##                      plot=T, ...) {
##  
##  if(stack_width > 1 | stack_width <= 0) stop('"stack_width" must be between 0 and 1')
##  if(is.null(data)) stop('Missing "data" parameter!')
##  
##  # Validate method
##  all_methods = c('bits', 'probability','custom')#, 'tsl')
##  pind = pmatch(method, all_methods)
##  method = all_methods[pind]
##  if(is.na(method)) stop("method must be one of 'bits' or 'probability', or 'custom'")
##  
##  # Convert character seqs to list
##  if(is.character(data) | is.matrix(data)) data = list("1"=data)
##  
##  if(is.list(data)){
##    # Set names for list if they dont exist
##    if(is.null(names(data))) names(data) = seq_along(data)
##    
##    lvls = names(data)
##    
##    # We have list of sequences - loop and rbind
##    data_sp = lapply(names(data), function(n){
##      curr_seqs = data[[n]]
##      logo_data(seqs = curr_seqs, method = method, stack_width = stack_width, 
##                rev_stack_order = rev_stack_order, seq_group = n, seq_type = seq_type, 
##                font = font, namespace=namespace)
##    })
##    data = do.call(rbind, data_sp)
##    # Set factor for order of facet
##    data$seq_group = factor(data$seq_group, levels = lvls)
##  }
##  
##  if(!plot) return(data)
##  
##  # Get sequence type
##  seq_type = attr(data, 'seq_type')
##  cs = get_col_scheme( col_scheme, seq_type )
##  
##  legend_title = attr(cs, 'cs_label')
##  
##  data = merge(data, cs, by='letter', all.x=T)
##  
##  # Make sure you retain order after merge
##  data = data[order(data$order),]
##  
##  # Do we have a gradient colscale
##  colscale_gradient = is.numeric( cs$group )
##  
##  colscale_opts = NULL
##  if(colscale_gradient){
##    # Set gradient colours 
##    colscale_opts = scale_fill_gradient(name=legend_title, low = low_col, 
##                                        high = high_col, na.value = na_col)
##  }else{
##    # Make group -> colour map
##    tmp = cs[!duplicated(cs$group) & !is.na(cs$group),]
##    col_map = unlist( split(tmp$col, tmp$group) )
##    
##    # Set colour scale options
##    colscale_opts = scale_fill_manual(values=col_map, name=legend_title, na.value=na_col)
##  } 
##  
##  # If letters and group are the same, don't draw legend
##  guides_opts = NULL
##  if(identical(cs$letter, cs$group)) guides_opts = guides(fill=F)
##  
##  y_lim = NULL
##  extra_opts = NULL
##  if(method == 'tsl'){
##    y_lab = 'Depleted    Enriched'
##    tmp = max(abs(data$y))
##    #y_lim = c(-tmp, tmp)
##    row_a = row_b = data[1,]
##    row_a$y = -tmp
##    row_b$y = tmp
##    data = rbind(data, row_a, row_b)
##    data$facet = factor(data$y > 0, c(T, F), c('Enriched', 'Depleted'))
##    extra_opts = NULL#facet_grid(facet~., scales='free')
##  }else if(method == 'custom'){
##    y_lab = ''
##  }else{
##    y_lab = method
##    substr(y_lab, 1, 1) = toupper(substr(y_lab, 1, 1))
##  }
##  
##  # Group data
##  data$group_by = with(data, interaction(seq_group, letter, position))
##  
##  data$x = data$x 
##  # Create layer
##  logo_layer = layer(
##    stat = 'identity', data = data, 
##    mapping = aes_string(x = 'x', y = 'y', fill='group', group='group_by'), 
##    geom = 'polygon', 
##    position = 'identity', show.legend = NA, inherit.aes = F,
##    params = list(na.rm = T, ...)
##  ) 
##  
##  
##  breaks_fun = function(lim){
##    # account for multiplicatuce expansion factor of 0.05
##    1: floor( lim[2] / 1.05 )
##  }
##  
##  # Expand 0.05 addidtive 
##  list(logo_layer, scale_x_continuous(breaks = breaks_fun, labels = identity), 
##       ylab(y_lab), xlab(''), colscale_opts, guides_opts, coord_cartesian(ylim=y_lim), 
##       extra_opts)
##}
##
##
###' Quick sequence logo plot
###' 
###' @description \code{ggseqlogo} is a shortcut for generating sequence logos. 
###' It adds the ggseqlogo theme \code{\link{theme_logo}} by default, and facets when multiple input data are provided. 
###' It serves as a convenient wrapper, so to customise logos beyond the defaults here, please use \code{\link{geom_logo}}.
###' 
###' @param data Character vector of sequences or named list of sequences. All sequences must have same width
###' @param facet Facet type, can be 'wrap' or 'grid'
###' @param scales Facet scales, see \code{\link{facet_wrap}}
###' @param ncol Number of columns, works only when \code{facet='wrap'}, see \code{\link{facet_wrap}}
###' @param nrow Number of rows, same as \code{ncol}
###' @param ... Additional arguments passed to \code{\link{geom_logo}}
###' 
###' @export
###' @examples
###' # Load sample data
###' data(ggseqlogo_sample)
###' 
###' # Plot a single DNA sequence logo
###' p1 = ggseqlogo( seqs_dna[[1]] )
###' print(p1)
###' 
###' # Plot multiple sequence logos at once
###' p2 = ggseqlogo( seqs_dna )
###' print(p2)
##ggseqlogo <- function(data, facet='wrap', scales='free_x', ncol=NULL, nrow=NULL, ...){
##  
##  # Generate the plot with default theme
##  p = ggplot() + geom_logo(data = data, ...) + theme_logo() 
##  
##  # If it's an inidivdual sequence logo, return plot
##  if(!'list' %in% class(data)) return(p)
##  
##  # If we have more than one plot, facet
##  facet_opts = c('grid', 'wrap')
##  pind = pmatch(facet, facet_opts)
##  facet = facet_opts[pind]
##  if(is.na(facet)) stop("facet option must be set to 'wrap' or 'grid'")
##  
##  if(facet == 'grid'){
##    p = p + facet_grid(~seq_group, scales = scales)
##  }else if(facet == 'wrap'){
##    p = p + facet_wrap(~seq_group, scales = scales, nrow = nrow, ncol = ncol)
##  }
##  
##  # Return plot
##  return(p)
##}
##
##
###' List of aligned transcription factor binding sequences 
###'
###' @name seqs_dna
###' @docType data
###' @keywords data
##NULL
##
###' List of aligned kinase-substrate binding sequences 
###'
###' @name seqs_aa
###' @docType data
###' @keywords data
##NULL
##
###' List of position frequency matrices for transcription factors
###'
###' @name pfms_dna
###' @docType data
###' @keywords data
##NULL
##
##
### message('-- running example')
### load('data/ggseqlogo_sample.rda')
### p = ggseqlogo(sample_data$seqs_dna, nrow=3)
### d = p$layers[[1]]$data
### print(p)
##
###' List color schemes available in ggseqlogo
###' 
###' @param v If true, font names are printed to stderr. Otherwise, color scheme names are returned as a character vector
###' @export
##list_col_schemes <- function(v=T){
##  
##  col_schemes = c('auto', 'chemistry', 'chemistry2','hydrophobicity', 'nucleotide', 'nucleotide2',
##             'base_pairing', 'clustalx', 'taylor')
##  if(!v) return(col_schemes)
##  message('Available ggseqlogo color schemes:')
##  for(f in col_schemes) message('\t', f)
##}
##
##
### Get color scheme
### @param col_scheme name of color scheme
### @param seq_type sequence type of color scheme
##get_col_scheme = function(col_scheme, seq_type='auto'){
##  
##  # Check if user-defined color scheme
##  if(is.data.frame(col_scheme)){
##    if(!'ggseqlogo_cs' %in% class(col_scheme)) 
##      stop('Colour scheme must be generated using "make_col_scheme" function')
##    return(col_scheme)
##  }
## 
##  # Get ambigious colour scheme
##  col_scheme = match.arg(col_scheme, list_col_schemes(F))
##  
##  # Get default color scheme for sequence type
##  if(col_scheme == 'auto'){
##    if(seq_type == 'auto') stop('"col_scheme" and "seq_type" cannot both be "auto"')
##    
##    col_scheme = switch(tolower(seq_type), aa = 'chemistry', 
##                        dna = 'nucleotide', rna = 'nucleotide', 
##                        other='nucleotide')
##
##  }
##  
##  
##  # Pick from default color schemes
##  cs = switch(col_scheme, 
##         # Color scheme based on chemistry of amino acids
##         chemistry2 = data.frame(
##           letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
##           group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
##           col = c(rep('#058644', 5), rep('#720091', 2), rep('#0046C5', 3), rep('#C5003E', 2), rep('#2E2E2E', 8)),
##           stringsAsFactors = F
##         ), 
##         
##         # Color scheme based on chemistry of amino acids
##         chemistry = data.frame(
##           letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
##           group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
##           col = c(rep('#109648', 5), rep('#5E239D', 2), rep('#255C99', 3), rep('#D62839', 2), rep('#221E22', 8)),
##           stringsAsFactors = F
##         ), 
##         
##         # Hydrophobicity index (PMID: 7108955) from -4.5 to 4.5
##         hydrophobicity = data.frame(
##           letter = c('I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'W', 
##                      'S', 'Y', 'P', 'H', 'D', 'E', 'N', 'Q', 'K', 'R'),
##           group = c(4.5, 4.2, 3.8, 2.8, 2.5, 1.9, 1.8, -0.4, -0.7, -0.9, -0.8,
##                       -1.3, -1.6, -3.2, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5),
##           stringsAsFactors=F
##         ), 
##         
##         # Colour based on nucleotide
##         nucleotide2 = data.frame(
##           letter = c('A', 'C', 'G', 'T', 'U'),
##           col = c('darkgreen', 'blue', 'orange', 'red', 'red'),
##           stringsAsFactors = F
##         ), 
##         
##         #alt red BA1200
##         nucleotide = data.frame(
##           letter = c('A', 'C', 'G', 'T', 'U'),
##           col = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839'),
##           stringsAsFactors = F
##         ), 
##         
##         base_pairing = data.frame(
##           letter = c('A', 'T', 'U', 'G', 'C'),
##           group = c(rep('Weak bonds', 3), rep('Strong bonds', 2)),
##           col = c(rep('darkorange', 3), rep('blue', 2)),
##           stringsAsFactors = F
##         ),
##         
##         # ClustalX color scheme: 
##         # http://www.jalview.org/help/html/colourSchemes/clustal.html
##         clustalx = data.frame(
##           letter = c('W', 'L', 'V', 'I', 'M', 'F', 'A', 'R', 'K', 'T', 'S', 'N', 'Q', 'D', 'E', 'H', 'Y', 'C', 'G', 'P'),
##           col = c(rep('#197FE5', 7), rep('#E53319', 2), rep('#19CC19', 4), rep('#CC4CCC', 2), 
##                   rep('#19B2B2', 2), '#E57F7F', '#E5994C', '#B0B000'),
##           stringsAsFactors = F
##         ),
##         
##         # Taylor color scheme (PMID: 9342138)
##         taylor = data.frame(
##           letter = c('D','S','T','G','P','C','A','V','I','L','M','F','Y','W','H','R','K','N','Q','E'),
##           col = c('#FF0000','#FF3300','#FF6600','#FF9900','#FFCC00','#FFFF00','#CCFF00','#99FF00',
##                   '#66FF00','#33FF00','#00FF00','#00FF66','#00FFCC','#00CCFF','#0066FF','#0000FF',
##                   '#6600FF','#CC00FF','#FF00CC','#FF0066'),
##           stringsAsFactors = F
##         )
##  )
##  
##  if(!'group' %in% names(cs)) cs$group = cs$letter
##  
##  # Set attributes
##  attr(cs, 'cs_label') = col_scheme
##  class(cs) = c('data.frame','ggseqlogo_cs')
##  
##  return(cs)
##}
##
##
##
##
##
###' Create new sequence logo color scheme
###' 
###' @param chars Vector of one letter characters 
###' @param groups Vector of groups for letters with same length as chars (optional if cols parameter is provided) 
###' @param cols Vector of colors with same length as chars (optional if values parameter is provided) 
###' @param values Vector of numerical values with same length as chars
###' @param name Name of color scheme
###' 
###' @export
###' 
###' @importFrom grDevices col2rgb
###' @examples 
###' 
###' # Discrete color scheme examples
###' cs1 = make_col_scheme(chars=c('A', 'T', 'G', 'C'), groups=c('g1', 'g1', 'g2', 'g2'), 
###'                       cols=c('red', 'red', 'blue', 'blue'), name='custom1')
###' 
###' cs2 = make_col_scheme(chars=c('A', 'T', 'G', 'C'), cols=c('red', 'red', 'blue', 'blue'), 
###'                       name='custom2')
###' 
###' # Quantitative color scheme
###' cs3 = make_col_scheme(chars=c('A', 'T', 'G', 'C'), values=1:4, name='custom3')
##make_col_scheme <- function(chars=NULL, groups=NULL, cols=NULL, values=NULL, name=''){
##  
##  
##  if(is.null(chars) | any(nchar(chars) != 1) | !is.character(chars))
##    stop('"chars" must be a character vector of one letter characters')
##  
##  
##  if(is.null(values)){
##    # Discrete colour scheme
##    
##    # Error check lengths
##    if(length(chars) != length(cols)) stop('"chars" and "cols" must have same length')
##    # Error check types
##    if(!is.character(cols)) stop('"cols" must be a character vector')
##    
##    # Check valid colours
##    tmp = col2rgb(cols); rm(tmp)
##    
##    if(is.null(groups)) groups = chars
##    
##    cs = data.frame( letter=chars, group=groups, col=cols, stringsAsFactors = F )
##    
##  }else{
##    
##    # Quantitative color scheme
##    if(length(chars) != length(values)) stop('"chars" and "values" must have same length')
##    cs = data.frame( letter=chars, group=values, stringsAsFactors=F )
##  }
##  
##  # Remove duplicate letters
##  cs = cs[!duplicated(cs$letter),]
##  
##  # Set attributes
##  attr(cs, 'cs_label') = name
##  class(cs) = c('data.frame','ggseqlogo_cs')
##  
##  return(cs)
##}


###############################################
#
# RWEBLOGO FUNCTIONS
#
###############################################



####################
#
# MOTIFX LAUNCH
#
##################

motifx <- function(fg.seqs, bg.seqs, central.res='ST', min.seqs=20, pval.cutoff=1e-6, verbose=T, perl.impl=F){
  
  AA = .GetAA()
  # Degenerate option not implemented yet.
  .CheckEmptySeqs <- function(){
    if(length(fg.seqs) == 0) stop('Could not find any foreground sequences!')
    if(length(bg.seqs) == 0) stop('Could not find any background sequences!')
  }
  
  .CheckEmptySeqs()
  
  cent.regex = ifelse(nchar(central.res) > 1, paste0('[',central.res,']'), central.res)
  c.res = strsplit(central.res, '')[[1]]
  c.res = intersect(AA, c.res)
  
  if(length(c.res) == 0) stop('Central residue must contain at least one amino acid character')
  
  # Check sequence widths 
  width = nchar(fg.seqs[1])
  if(width < 3 | width > 35) stop('Sequence width must be between 3 and 35!')
  if(width %% 2 == 0) stop('Sequence width must be an odd number!')
  if(width != nchar(bg.seqs[1])) 
    stop('Widths for foreground and background data must be equal!')
  
  # Ensure k-mers have same lengths
  nc.pos = nchar(fg.seqs)
  if( any(nc.pos[1] != nc.pos[-1]) ) stop('Foreground k-mers must be same lentth.')
  nc.bg = nchar(bg.seqs)
  if( any(nc.bg[1] != nc.bg[-1]) ) stop('Background k-mers must be same lentth.')
  
  # Get central index
  ci = ceiling(width/2)
  
  # Filter the central residue for allowed residues
  fg.seqs = fg.seqs[substr(fg.seqs, ci,ci) %in% c.res]
  bg.seqs = bg.seqs[substr(bg.seqs, ci,ci) %in% c.res]
  
  
  if(T){
    # Only set to true for exact match to motifx webserver 
    # Remove non-amino acid residues
    fg.seqs = fg.seqs[!grepl('\\-|\\*|[BJOUXZ]', fg.seqs)]
    bg.seqs = bg.seqs[!grepl('\\-|\\*|[BJOUXZ]', bg.seqs)]
    bg.seqs = unique(bg.seqs)
  }
  
  .CheckEmptySeqs()
  
  seqs = fg.seqs
  data = list()
  while(TRUE){
    # Find the motif
    mt = .FindMotif(fg.seqs = seqs, bg.seqs = bg.seqs, min.seqs = min.seqs, 
                   pval.cutoff = pval.cutoff, cent.regex = cent.regex, verbose = verbose, 
                   perl.impl=perl.impl)
    if(is.null(mt)) break

		####STEF create file with sequence specific for motif
		write.table(mt$pos,file=paste("sequenceMotif_",mt$motif.regex,"_.txt",sep=""), col.names = F, row.names = F)
	
    # Append to list of data
    data = append(data, list(mt))
    # Remove stuff already apart of a motif
    if(verbose) writeLines('Step 2: positive and negative set reduction')
    seqs = seqs[! grepl(mt$motif.regex, seqs)]
    bg.seqs = bg.seqs[! grepl(mt$motif.regex, bg.seqs)]
    
    # No more sequences left to process, breaking.
    if(length(seqs) < min.seqs) break

		##STEF create Logo
		##weblogo(seqs=mt$pos, format="png", file.out=paste("logoMotif_",mt$motif.regex,"_.png",sep=""),verbose=T,size='large', show.xaxis = F, show.yaxis = F, errorbars = F, reverse.stacks =T,color.scheme = 'chemistry')
		
		##create LOGO woth ggseqlogo
		png(filename=paste("logoMotif_",mt$motif.regex,"_.png",sep=""),width=800, height = 400, unit="px")
		plot(ggseqlogo(mt$pos,seq_type='aa'))
		dev.off()
		
	
	}
  
  if(verbose) writeLines('Converged, no more enrichments!')
  
  df = data.frame(motif = sapply(data, function(x) x$motif.regex), 
                  score = sapply(data, function(x) x$motif.score), 
                  fg.matches = sapply(data, function(x) x$fg.matches), 
                  fg.size = sapply(data, function(x) x$fg.size), 
                  bg.matches = sapply(data, function(x) x$bg.matches), 
                  bg.size = sapply(data, function(x) x$bg.size),
                  stringsAsFactors=F)
  
  df$fold.increase = (df$fg.matches/df$fg.size)/(df$bg.matches/df$bg.size)
  rownames(df) = NULL
  if(nrow(df) == 0) return(NULL)
  return(df)
}

####>Revision history<####
# 1.0.0 first version (SL 26/07/17)
