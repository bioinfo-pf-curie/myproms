################################################################################
# correlation.R      1.0.0                                                     #
# Authors: G. Arras,P. Poullet (Institut Curie)                                #
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

# Usage:
# R CMD BATCH --no-save --no-restore "--args <correlation/distance> <fileName>" correlation.R

correlation <- function(filename) {
  data <- read.table(file=filename,sep="\t",header=T,row.names=1)

  cor.pearson <- cor(data,use="pairwise.complete.obs",method="pearson")
  cor.spearman <- cor(data,use="pairwise.complete.obs",method="spearman")

  write.table(cor.pearson,file=paste("Results_Pearson_",filename,sep=""),row.names=F,quote=F,sep="\t")
  write.table(cor.spearman,file=paste("Results_Spearman_",filename,sep=""),row.names=F,quote=F,sep="\t")
}

distance <- function(filename) {
  data <- read.table(file=filename,sep="\t",header=T,row.names=1)
  data <- t(data)

  #dist.manhattan <- dist(data,method="manhattan",diag=T,upper=T)
  dist.binary <- dist(data,method="binary",diag=T,upper=T)

  #write.table(as.matrix(dist.manhattan),file=paste("Results_Manhattan_",filename,sep=""),row.names=F,quote=F,sep="\t")
  write.table(as.matrix(dist.binary),file=paste("Results_Binary_",filename,sep=""),row.names=F,quote=F,sep="\t")
}

#> Main
args <- commandArgs(trailingOnly=TRUE)

method <- as.character(args[1])
dataFile <- as.character(args[2])

if (method=='correlation') {
  correlation(dataFile)
} else {
  distance(dataFile)
}

####> Revision history
# 1.0.0 First stable version (PP 06/02/18)
