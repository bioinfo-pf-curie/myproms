################################################################################
# prepareExplorAna.R       1.0.6                                               #
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
# Impute +/- infinite ratio, missing values using the missMDA package (PCA-base imputation) & filters rows based on top N larger sd

#English
Sys.setenv( LANGUAGE = "en")

paramR<-read.table("R_parameters.txt", header=T, sep="\t")
inputMat <- read.table("matrix.txt", header=TRUE, check.names=FALSE)

impMat=data.frame()

#write.table(inputMat,file="matrixWithMissingValues.txt",quote=FALSE, sep="\t", col.names=FALSE)

#####----- Replacement of +/- infinite ratios ------#####
# If exist infinite ratios (+/-1000 values) then they are NOT treated as NA (otherwise already replaced by NA by startExploratoryAnalysis.cgi)
### Compute st.dev, mean, ref max, ref min
mat.sd <- sd(inputMat[!is.na(inputMat) & (inputMat != 1000) & (inputMat != -1000)])
mat.mean <-mean(inputMat[!is.na(inputMat) & (inputMat!= 1000)& (inputMat != -1000)])
mat.min <- mat.mean - 4 * mat.sd
mat.max <- mat.mean + 4 * mat.sd

### Replace +/-inf values with random values in the range of mean +/- 4 st.dev +/- random(+/- 0.5 st.dev)
inputMat <- apply(inputMat,1, function(x) { replace(x,which(x == 1000),runif(length(which(x== 1000)), mat.max-0.5*mat.sd,mat.max+0.5*mat.sd)) })
inputMat <- t(inputMat)
inputMat <- apply(inputMat,1, function(x) { replace(x,which(x == -1000),runif(length(which(x== -1000)), mat.min-0.5*mat.sd,mat.min+0.5*mat.sd)) })
inputMat <- t(inputMat)


#####----- NA imputation using missMDA ------#####
#numbNA <- length(inputMat[is.na(inputMat)])
numbNA <- length(which(is.na(inputMat)))

##Launch missMDA if missing values
if (paramR$MISSING_VALUE != 0 && numbNA > 0) { # missMDA fails if no NA at all

	library(missMDA)
	### Estimate number of composantes required
	nb <- estim_ncpPCA(inputMat,ncp.max=10)
	nb$ncp

	### Replacement of NAs
	resImpute <- imputePCA(inputMat,ncp=nb$ncp)

	impMat <- resImpute$completeObs
    
    ## Create distribution plot for existing/missing values
    library(plyr)
    library(ggplot2)
    theme_set(theme_classic())
    
    # Categorize missing and imputed values
    rep = rbind(data.frame(type=factor('Existing values'), value=impMat[!is.na(inputMat)]), data.frame(type=factor('Imputed'), value=impMat[is.na(inputMat)]))
    
    # Get mean of each set 
    mu <- ddply(rep, "type", summarise, grp.mean=mean(value))
    
    # Plot both density distribution
    p <- ggplot(rep, aes(value, fill=type, color=type)) +
      geom_density(aes(y = ..count..), alpha=.2) +
      geom_vline(data=mu, aes(xintercept=grp.mean, color=type), 
                 linetype="dashed", size=0.5) + 
      labs(title="Values distribution", 
           subtitle="For both initial and imputed ones",
           x="Quantification value",
           y="Amount of proteins")
    
    png(filename="valueDistribution.png")
    plot(p)
    dev.off()
} else {
	impMat<-inputMat
}

##Modification pipeline, exclude ambiguity isoforms
if (paramR$EXCLUDE_AMB == "TRUE") {
    impMat<-impMat[grep(":",rownames(impMat),invert=T),]
}

if (paramR$PROTEIN_SELECTION == "none") {
	write.table(impMat,file="matrixProcessed.txt",quote=FALSE, sep="\t", col.names=NA)
} else {
	prot.sd=matrix()
    for (i in 1:length(impMat[,1])) {
		  prot.sd[i]=sd(impMat[i,],na.rm=TRUE)
    }

    prot.sd<-as.data.frame(prot.sd)
    rownames(prot.sd)<-rownames(impMat)
	impMat_merge<-impMat[order(prot.sd$prot.sd, decreasing=TRUE),]
    nbRowSD<-nrow(impMat_merge)

    if (paramR$AGGREGATE == "TRUE") {
  		prot.sd$PROTEIN<-gsub("-.*","", rownames(prot.sd))
  		impMat_merge<-impMat_merge[!duplicated(prot.sd$PROTEIN),]
  		nbRowSD<-nrow(impMat_merge)
    }

    if (paramR$KEEP_PROT == "TRUE") {
		  nbProtSD<-ifelse(paramR$KEEP_PROT_NB > nbRowSD, nbRowSD, paramR$KEEP_PROT_NB)
		  impMat_merge<-impMat_merge[1:nbProtSD,]
    }

	write.table(impMat_merge,file="matrixProcessed.txt",quote=FALSE, sep="\t", col.names=NA)
}

####>Revision history<####
# 1.0.6 [BUGFIX] NA distribution plot is not drawn is no NA (PP 13/11/19)
# 1.0.5 [FEATURE] Add density plot (count normalized) as imputation quality representation (VS 14/10/19)
# 1.0.4 rename AGREGATE to AGGREGATE and comment ambiguity exclusion, manage by perl(SL 28/09/18)
# 1.0.3 add modification pipeline (SL 28/09/17)
# 1.0.2 Moved nb$ncp within definition scope of nb (PP 11/01/17)
# 1.0.1 TopN proteins based on sd & +/-infinite imputations (PP 18/02/16)
# 1.0.0 Former imputeNA.R (PP 21/01/16)
