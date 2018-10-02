################################################################################
# goQuanticluster.R       1.0.0                                                #
# Authors: Patrick Poullet, Stephane Liva, Guillaume Arras (Institut Curie)    #
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
#library("gplots")
args <- commandArgs(TRUE)
mergeFile <- args[1]
data <- read.table(mergeFile,header=F,row.names=1,sep="\t")
data.log <- log10(as.matrix(data))
colnames(data.log) <- c("bin1","bin2","bin3","bin4","bin5")
# Z-scoring
for(i in 1:nrow(data.log)){
   mean.l <- mean(data.log[i,])
   sd.l <- sd(data.log[i,])
   data.log[i,] <- (data.log[i,]-mean.l)/sd.l
}
clust.prot <- hclust(dist(data.log), method="average")
write.table(data.frame(clust.prot$merge,sort(clust.prot$height)),file="goDendro.txt",quote=FALSE, sep="\t", col.names=NA )
write.table(data.frame(clust.prot$order,clust.prot$labels[clust.prot$order]), file="goOrder.txt",quote=FALSE, sep="\t", col.names=NA )

#hmCarpet <- heatmap.2(as.matrix(data.log), distfun = function(x) dist(x,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'average'), key=T, symkey=FALSE, density.info="none", trace="none")$carpet
#writeLines(colnames(hmCarpet), con="extractOrder.txt")


####>Revision history<####
# 1.0.0 First version, replace R cluster from displayGOQuantiAnalysis, replace heatmap.2 by hclust (SL 05/12/14)
