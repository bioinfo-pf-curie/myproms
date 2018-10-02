################################################################################
# cluster.R       2.0.3                                                        #
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

#English
Sys.setenv( LANGUAGE = "en")

args <- commandArgs(TRUE)
cluster.metric <- args[1]
cluster.method <- args[2]
metricType <-args[3]

Matrix <- read.table("matrix.txt", header=TRUE, check.names=FALSE)

MatrixTranspo <- t(Matrix)
# quantifications: transpo, proteins: non transpo
if (metricType == "cor") {
    metric.protein <- as.dist(1-cor(MatrixTranspo, method=cluster.metric))
    metric.quantif <- as.dist(1-cor(Matrix, method=cluster.metric))

    #c.protein<-hclust(as.dist(1-cor(MatrixTranspo)), method=cluster.method)
    #c.quantif<-hclust(as.dist(1-cor(Matrix)), method=cluster.method)
}else {
    metric.protein <- dist(Matrix, method=cluster.metric)
    metric.quantif <- dist(MatrixTranspo, method=cluster.metric)
    #c.protein<-hclust(as.dist(1-cor(MatrixTranspo)), method="centroid")
    #c.quantif<-hclust(as.dist(1-cor(Matrix)), method="centroid")
}
c.protein<-hclust(metric.protein, method=cluster.method)
c.quantif<-hclust(metric.quantif, method=cluster.method)
#c.protein<-clustering(MatrixTranspo, metric=cluster.metric, method=cluster.method)
#c.quantif<-clustering(as.matrix(Matrix), metric=cluster.metric, method=cluster.method)

write.table(data.frame(c.protein$order,c.protein$labels[c.protein$order]), file="protCluster.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(data.frame(c.protein$merge,sort(c.protein$height)), file="protDendro.txt",quote=FALSE, sep="\t", col.names=NA )
write.table(data.frame(c.quantif$order,c.quantif$labels[c.quantif$order]), file="quantifCluster.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(data.frame(c.quantif$merge,sort(c.quantif$height)), file="quantifDendro.txt",quote=FALSE, sep="\t", col.names=NA )

write.table(MatrixTranspo, file="MatrixTranspo.txt", quote=FALSE, sep="\t", col.names=NA)

####>Revision history<####
# 2.0.3 +/-infinite imputations moved to prepareExplorAna.R (PP 18/02/16)
# 2.0.2 Bug fix in mysd calculation (PP 03/02/16)
# 2.0.1 Replace +/- infinite ratios [-1000,1000] & sets language to English (PP 21/01/16) (PP 29/10/15)
# 2.0.0 new version of clustering using Hclust, detect correlation or distance metric (SL 17/12/14)
# 1.0.1 Rename R objects and files (PP 07/08/14)
# 1.0.0 First version (GA 08/07/14)
