################################################################################
# PCA.R       1.0.8                                                            #
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

require(FactoMineR)

Matrix <- read.table("matrix.txt", header=TRUE, check.names=FALSE)

MatrixTranspo <- t(Matrix)
write.table(MatrixTranspo,file="MatrixTranspo.txt",quote=FALSE, sep="\t", col.names=NA)

ProtPCASc <- PCA(as.data.frame(Matrix),scale.unit=TRUE, graph=FALSE, ncp=min(dim(Matrix)))
QuantifPCASc <- PCA(as.data.frame(MatrixTranspo),scale.unit=TRUE, graph=FALSE, ncp=min(dim(Matrix)))
ProtPCA <- PCA(as.data.frame(Matrix),scale.unit=FALSE, graph=FALSE, ncp=min(dim(Matrix)))
QuantifPCA <- PCA(as.data.frame(MatrixTranspo),scale.unit=FALSE, graph=FALSE, ncp=min(dim(Matrix)))

# Scaled
write.table(ProtPCASc$ind$coord,file="protCoordinates_sc.txt", quote=FALSE, sep="\t", col.names=NA)

ProtPCASc$eig<-as.data.frame(ProtPCASc$eig) #to be compatible with previous version of FactoMineR
write.table(ProtPCASc$eig$"percentage of variance",file="protContrib_sc.txt", quote=FALSE, sep="\t", col.names=NA)

# Scaled transpose
write.table(QuantifPCASc$ind$coord,file="quantifCoordinates_sc.txt", quote=FALSE, sep="\t", col.names=NA)

QuantifPCASc$eig<-as.data.frame(QuantifPCASc$eig) #to be compatible with previous version of FactoMineR
write.table(QuantifPCASc$eig$"percentage of variance",file="quantifContrib_sc.txt", quote=FALSE, sep="\t", col.names=NA)

# Not scaled
write.table(ProtPCA$ind$coord,file="protCoordinates.txt", quote=FALSE, sep="\t", col.names=NA)

ProtPCA$eig<-as.data.frame(ProtPCA$eig) #to be compatible with previous version of FactoMineR
write.table(ProtPCA$eig$"percentage of variance",file="protContrib.txt", quote=FALSE, sep="\t", col.names=NA)

# Not scaled transpose
write.table(QuantifPCA$ind$coord,file="quantifCoordinates.txt", quote=FALSE, sep="\t", col.names=NA)

QuantifPCA$eig<-as.data.frame(QuantifPCA$eig) #to be compatible with previous version of FactoMineR
write.table(QuantifPCA$eig$"percentage of variance",file="quantifContrib.txt", quote=FALSE, sep="\t", col.names=NA)

# Dimension description
nb.axes <- min(dim(Matrix))-1
QuantifPCASc.dimdesc <- dimdesc(QuantifPCASc, axes=1:nb.axes, proba = 0.05) #proba: the significance threshold considered to characterized the dimension
for(i in 1:nb.axes)
{
#    write.table(QuantifPCASc.dimdesc[[paste("Dim",i,sep=".")]],file=paste(paste("quantifProtDim",i,"_sc",sep=""),"txt",sep="."), quote=FALSE, sep="\t", col.names=NA)
    write.table(na.omit(QuantifPCASc.dimdesc[[paste("Dim",i,sep=".")]][["quanti"]]),file=paste(paste("quantifProtDim",i,"_sc",sep=""),"txt",sep="."), quote=FALSE, sep="\t", col.names=NA)
}

QuantifPCA.dimdesc <- dimdesc(QuantifPCA, axes=1:nb.axes, proba = 0.05)
for(i in 1:nb.axes)
{
#    write.table(QuantifPCA.dimdesc[[paste("Dim",i,sep=".")]],file=paste(paste("quantifProtDim",i,sep=""),"txt",sep="."), quote=FALSE, sep="\t", col.names=NA)
    write.table(na.omit(QuantifPCA.dimdesc[[paste("Dim",i,sep=".")]][["quanti"]]),file=paste(paste("quantifProtDim",i,sep=""),"txt",sep="."), quote=FALSE, sep="\t", col.names=NA)
}

####>Revision history<####
# 1.0.8 compatible with previous version of factoMineR (SL 18/10/17)
# 1.0.7 +/-infinite imputations moved to prepareExplorAna.R (PP 18/02/16)
# 1.0.6 Bug fix in mysd calculation (PP 03/02/16)
# 1.0.5 Sets language to English (PP 21/01/16)
# 1.0.4 Prevents NAs and row index to be printed to quantifProtDim files (PP 10/11/15)
# 1.0.3 Replace +/- infinite ratios [-1000,1000] (PP 29/10/15)
# 1.0.2 Rename R objects and files (PP 07/08/14)
# 1.0.1 Change file names (GA 08/07/14)
# 1.0.0 First version (GA 08/07/14)
