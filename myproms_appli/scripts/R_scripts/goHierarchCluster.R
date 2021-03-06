################################################################################
# goHierarchCluster.R       1.0.0                                              #
# Authors: Valentin Sabatet (Institut Curie)                                   #
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
args <- commandArgs(TRUE)
binMatrixFile <- args[1]
maxGOTerms <- args[2]
maxProts <- args[3]
data <- read.table(binMatrixFile,header=T,row.names=1,sep="\t")

d <- dist(data, method = "binary")
clust.prot <- hclust(d, method="ward.D")
write.table(data.frame(clust.prot$merge,sort(clust.prot$height)),file= paste("goDendro_", maxGOTerms, "_", maxProts, ".txt", sep="") ,quote=FALSE, sep="\t", col.names=NA )
write.table(data.frame(clust.prot$order,clust.prot$labels[clust.prot$order]), file= paste("goOrder_", maxGOTerms, "_", maxProts, ".txt", sep=""),quote=FALSE, sep="\t", col.names=NA )

####>Revision history<####
# 1.0.0 First version (VS 08/01/20)
