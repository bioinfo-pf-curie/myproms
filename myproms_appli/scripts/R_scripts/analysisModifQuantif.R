################################################################################
# analysisModifQuantif.R       1.1.0                                           #
# Authors: Patrick Poullet, Stephane Liva, Pierre Gestraud (Institut Curie)    #
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

##delete previous session
rm(list=ls())
Sys.setenv( LANGUAGE = "en")

args <- commandArgs(TRUE)
matrixExport <- args[1]
pepFileExport <-args[2]

## Load R Param file from myProms
paramR<-read.table("R_parameters.txt", header=T, sep="\t")

#########Load matrix from myProMS
matrixPTM<-read.table(matrixExport, header=T, sep="\t",row.names=1, check.names=F)
##matInfoValue<-as.matrix(read.table("value_info.txt",header=T,sep="\t",row.names=1, check.names=F,stringsAsFactors = F))
##existSdGeo <- FALSE
##if (file.exists("sd_geo.txt")) {
##	matSDgeo<-read.table("sd_geo.txt",header=T,sep="\t",row.names=1, check.names=F,stringsAsFactors = F)
##	existSdGeo <- TRUE
##}

#####----- Impute missing values ------#####
impMat=data.frame()
if (paramR$IMPUTE == "TRUE") {

	#####----- Replacement of +/- infinite ratios (+/-1000 values)------#####
	#if (!is.na(as.vector(table(abs(matrixPTM)==1000))[2])) { #} FALSE TRUE matrix (only 1 type if 100% FALSE/TRUE but cannot be 100%==1000!)
	if (length(which(abs(matrixPTM) == 1000)) > 0) {	
		
		### Compute st.dev, mean, ref max, ref min
		mat.sd <- sd(matrixPTM[!is.na(matrixPTM) & (matrixPTM != 1000) & (matrixPTM != -1000)])
		mat.mean <-mean(matrixPTM[!is.na(matrixPTM) & (matrixPTM!= 1000) & (matrixPTM != -1000)])
		mat.min <- mat.mean - 4 * mat.sd
		mat.max <- mat.mean + 4 * mat.sd
		
		### Replace +/-inf values with random values in the range of mean +/- 4 st.dev +/- random(+/- 0.5 st.dev)
		matrixPTM <- apply(matrixPTM,1, function(x) { replace(x,which(x == 1000),runif(length(which(x== 1000)), mat.max-0.5*mat.sd,mat.max+0.5*mat.sd)) })
		matrixPTM <- t(matrixPTM)
		matrixPTM <- apply(matrixPTM,1, function(x) { replace(x,which(x == -1000),runif(length(which(x== -1000)), mat.min-0.5*mat.sd,mat.min+0.5*mat.sd)) })
		matrixPTM <- t(matrixPTM)
	}
	
	#####----- NA imputation using missMDA ------#####
	if (length(matrixPTM[is.na(matrixPTM)]) > 0) { # Fails if no NA
		library(missMDA)
		
		### Estimate number of composantes required
		nb <- estim_ncpPCA(matrixPTM,ncp.max=10)
		nb$ncp
	
		### Replacement of NAs
		resImpute <- imputePCA(matrixPTM,ncp=nb$ncp)
	
		impMat <- resImpute$completeObs
	} else {
		impMat <- matrixPTM
	}
} else {
	impMat<-matrixPTM
}

#####----- Filter data set ------#####
filteredData <- FALSE
PTMSign=data.frame()
##NONE
if (paramR$PROTEIN_SELECTION == "none") {
	PTMSign<-impMat
	##if (paramR$FOCUS != 0){
	##	PTMSign$SITE<-rownames(impMat)
	##} else {
	##	PTMSign$PROTEIN<-rownames(impMat)
	##}
}

##SAMPLE
if (paramR$PROTEIN_SELECTION == "sample") {
	filteredData <- TRUE
    prot.sd=matrix()
    for (i in 1:length(impMat[,1])) {
		prot.sd[i]=sd(impMat[i,],na.rm=TRUE)
    }

    prot.sd<-as.data.frame(prot.sd)
    rownames(prot.sd)<-rownames(impMat)
    colnames(prot.sd)<-"sd"
	if (paramR$FOCUS != 0){
		prot.sd$PROTEIN<-gsub("-.*","", rownames(prot.sd))
		prot.sd$SITE<-rownames(prot.sd)
	} else {
		prot.sd$PROTEIN<-rownames(prot.sd)
	}
    prot.sd<-na.omit(prot.sd)
    nbRowSD<-nrow(prot.sd)

    if (paramR$AGGREGATE == "TRUE") {
		 prot.sd<-prot.sd[order(prot.sd$sd, decreasing=TRUE),]
		 prot.sd<-prot.sd[!duplicated(prot.sd$PROTEIN),]
		 nbRowSD<-nrow(prot.sd)
    }

    if (paramR$KEEP_PROT == "TRUE") {
		if (paramR$KEEP_PROT_NB == -1) {
			nbProtSD<-nbRowSD
		} else {
			nbProtSD<-ifelse(paramR$KEEP_PROT_NB > nbRowSD, nbRowSD, paramR$KEEP_PROT_NB)
		}
		prot.sd<-prot.sd[1:nbProtSD,]
    }

    PTMSign<-prot.sd
	##if (paramR$FOCUS != 0){
	##	impMat<-impMat[PTMSign$SITE,,drop=FALSE]
	##} else {
	##	impMat<-impMat[PTMSign$PROTEIN,,drop=FALSE]
	##}
	
	### Clean up dataframe
	PTMSign$PROTEIN <- NULL
	PTMSign$SITE <- NULL
	PTMSign$sd <- NULL
}

##GROUP
if (paramR$PROTEIN_SELECTION == "group") {
	filteredData <- TRUE
	
    #############load group info
    group<-read.table("groups.txt",header=T,sep="\t",row.names=1)

    ##check if there is at least 2 values by group
    nbValuesByProtBygroup<-apply(impMat,1, function(x) aggregate(x, by=list(group=group$GROUP),FUN=function(y) sum(!is.na(y))>=2))

    ##get all values
    allValues<-sapply(nbValuesByProtBygroup, FUN = function(x) all(x$x) )

    ##get TRUE values
    trueValues<-as.matrix(which(allValues))

    impMat<-impMat[rownames(trueValues),]

    #launch ANOVA FOR GROUP
    matrixPTM_Anova<-as.data.frame(apply(impMat, 1,function(x) anova(lm(x ~ group$GROUP))[1,5]))
    colnames(matrixPTM_Anova)<-"pvalue"
    if (paramR$FOCUS != 0){
		matrixPTM_Anova$PROTEIN<-gsub("-.*","", rownames(matrixPTM_Anova))
		matrixPTM_Anova$SITE<-rownames(matrixPTM_Anova)
	} else {
		matrixPTM_Anova$PROTEIN<-rownames(matrixPTM_Anova)
	}
    matrixPTM_Anova<-na.omit(matrixPTM_Anova)
    nbRowAnova<-nrow(matrixPTM_Anova)

    ##P-value threshold for Anova
    if (paramR$ANOVA_P_VALUE_CHK == "TRUE") {
		matrixPTM_Anova<-matrixPTM_Anova[matrixPTM_Anova <= paramR$ANOVA_P_VALUE,,drop=FALSE]
		nbRowAnova<-nrow(matrixPTM_Anova)
    }

    ##Best sites for a protein
    if (paramR$AGGREGATE == "TRUE") {
		matrixPTM_Anova<-matrixPTM_Anova[order(matrixPTM_Anova$pvalue, decreasing=FALSE),]
		matrixPTM_Anova<-matrixPTM_Anova[!duplicated(matrixPTM_Anova$PROTEIN),,drop=FALSE]
		nbRowAnova<-nrow(matrixPTM_Anova)
    }

    ##Keep prot
    if (paramR$KEEP_PROT == "TRUE") {
		nbProtAnova<-ifelse(paramR$KEEP_PROT_NB > nbRowAnova, nbRowAnova, paramR$KEEP_PROT_NB)
		matrixPTM_Anova<-matrixPTM_Anova[1:nbProtAnova,,drop=FALSE]
    }

    PTMSign<-matrixPTM_Anova
	##if (paramR$FOCUS != 0){
	##	impMat<-impMat[PTMSign$SITE,,drop=FALSE]
	##} else {
	##	impMat<-impMat[PTMSign$PROTEIN,,drop=FALSE]
	##}

	### Clean up dataframe
	PTMSign$PROTEIN <- NULL
	PTMSign$SITE <- NULL
	PTMSign$pvalue <- NULL
}


#######load annotation
##matrixAnnot<-read.csv("annotation.txt",sep="\t", header=T, quote = "", stringsAsFactor=FALSE)
##
##if (paramR$FOCUS != 0){## PTM site
##	listAnnotPTM<-cbind(rownames(impMat), gsub("-.*","",rownames(impMat)))
##	colnames(listAnnotPTM)<-c("SITE","PROTEIN")
##	annotPTM<-merge(matrixAnnot, listAnnotPTM, by.x="PROTEIN", by.y="PROTEIN")
##	annotPTM$SITE<-as.character(annotPTM$SITE)
##} else {
##	listAnnotPTM<-as.data.frame(rownames(impMat))
##	colnames(listAnnotPTM)<-"PROTEIN"
##	annotPTM<-merge(matrixAnnot, listAnnotPTM, by.x="PROTEIN", by.y="PROTEIN")
##	annotPTM$PROTEIN<-as.character(annotPTM$PROTEIN)
##}
##existPvalue <- FALSE
##if (paramR$QUANTIF_FAM == "RATIO" && paramR$FOCUS != 2) {
##	#####load matrix pep
##	matrixPepPTM<-read.table(pepFileExport,sep="\t",header=T,row.names=1,check.names=F)
##	
##	#####load matrix pValue
##	if (file.exists("pvalue.txt")) {
##		existPvalue <- TRUE
##		matrixPvaluePTM<-read.table("pvalue.txt",sep="\t",header=T,row.names=1,check.names=F)
##
##		if (paramR$FOCUS != 0){
##			pepPTM<-matrixPepPTM[PTMSign$SITE,,drop=FALSE]
##			pValuePTM<-matrixPvaluePTM[PTMSign$SITE,,drop=FALSE]
##		} else {
##			pepPTM<-matrixPepPTM[PTMSign$PROTEIN,,drop=FALSE]
##			pValuePTM<-matrixPvaluePTM[PTMSign$PROTEIN,,drop=FALSE]
##		}
##	}
##}


write.table(PTMSign,paste("processed_",matrixExport,sep=""), col.names=NA, quote=F,sep="\t")

if (filteredData == TRUE) {
	extraFiles<-c('annotation.txt',pepFileExport,'pvalue.txt','sd_geo.txt','value_info.txt')
	for (exFile in extraFiles) {
		if (file.exists(exFile)) {
			if (exFile=='annotation.txt') {
				exData<-read.csv(exFile,sep="\t", header=T, quote = "", stringsAsFactor=FALSE)
			} else if (exFile=='value_info.txt') { # requires special handling
				exData<-as.matrix(read.table(exFile,header=T,sep="\t",row.names=1, check.names=F,stringsAsFactors = F))
			} else {
				exData<-read.table(exFile,sep="\t",header=T,row.names=1,check.names=F)
			}
			exData<-exData[rownames(PTMSign),]
			###if (paramR$FOCUS != 0) {
			###	exData<-exData[PTMSign$SITE,,drop=FALSE]
			###} else {
			###	exData<-matInfoValue[PTMSign$PROTEIN,,drop=FALSE]
			###}
			write.table(exData,paste("processed_",exFile,sep=""),col.names=NA, quote=F,sep="\t")
			
		}
	}
}

## Heatmap of data type
if (paramR$IMPUTE == "TRUE") {
	valueFile <- 'value_info.txt'
	if (filteredData == TRUE) {
		valueFile <- 'processed_value_info.txt'
	}
	if (file.exists(valueFile)) {
		infoValue<-as.matrix(read.table(valueFile,header=T,sep="\t",row.names=1, check.names=F,stringsAsFactors = F))
		infoValue[infoValue == "Val"]<- 2
		infoValue[infoValue == "MV"]<- 0
		infoValue[infoValue == "-Inf"]<- -1
		infoValue[infoValue == "+Inf"]<- 1
		
		write.table(infoValue,"processed_value_info_sub.txt",quote=F,sep="\t",col.name=NA)
		
		infoHM<-as.matrix(read.table("processed_value_info_sub.txt", header=T, row.names=1,check.names=1,check.names=F))
		
		require(pheatmap)
		png("processed_value_info.png", unit='in', width = 5, height = 15, res=600)
		pheatmap(infoHM,cluster_rows = F,cluster_cols = F, show_rownames = F, color=c("blue","black","red","white"),legend_breaks = c(-1,0,1,2),legend_labels = c("-Infinite","Missing values","+Infinite","Values"),border_color=NA)
		dev.off()
	}
}
#if (paramR$FOCUS != 0){
#	infoValue<-matInfoValue[PTMSign$SITE,,drop=FALSE]
#} else {
#	infoValue<-matInfoValue[PTMSign$PROTEIN,,drop=FALSE]
#}

##if (paramR$FOCUS != 0) {
##	annotPTM<-merge(annotPTM,PTMSign,by.x="SITE", by.y="SITE")
##} else {
##	annotPTM<-merge(annotPTM,PTMSign,by.x="PROTEIN", by.y="PROTEIN")
##}

##if (paramR$QUANTIF_FAM == "RATIO") {
##	if (paramR$FOCUS != 2){
##		write.table(pepPTM,paste("processed_",pepFileExport,sep=""),col.names=NA, quote=F,sep="\t")
##		if (existPvalue) {
##			write.table(pValuePTM,"processed_pvalue.txt",col.names=NA, quote=F,sep="\t")
##		}
##	}
##	write.table(impMat,"processed_ratio.txt", col.names=NA, quote=F,sep="\t")
##} else {
##	write.table(impMat,paste("processed_",matrixExport,sep=""), col.names=NA, quote=F,sep="\t")
##}

##if (paramR$FOCUS != 0){##modification
##	write.table(annotPTM[,c(1,2,3,4,5,6)],"processed_annotation.txt",sep="\t",quote=F, row.names = F)
##} else {
##	write.table(annotPTM[,c(1,2,3,4,5)],"processed_annotation.txt",sep="\t",quote=F, row.names = F)
##}

if (paramR$PROTEIN_SELECTION == "sample") {
	if (paramR$FOCUS != 0){##modification
		write.table(cbind(annotPTM$GENE_NAME,annotPTM$SITE,annotPTM$sd),"sd.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE_NAME","SITE","SD"))
	} else {
		write.table(cbind(annotPTM$GENE_NAME,annotPTM$PROTEIN,annotPTM$sd),"sd.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE_NAME","PROTEIN","SD"))
	}
} else if (paramR$PROTEIN_SELECTION == "group") {
	if (paramR$FOCUS != 0){##modification
		write.table(cbind(annotPTM$GENE_NAME,annotPTM$SITE,annotPTM$pvalue),"anova_pvalue.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE_NAME","SITE","PVALUE"))
	} else {
		write.table(cbind(annotPTM$GENE_NAME,annotPTM$PROTEIN,annotPTM$pvalue),"anova_pvalue.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE_NAME","PROTEIN","PVALUE"))
	}
}


##if (paramR$FOCUS != 0) {
##	infoValue<-matInfoValue[PTMSign$SITE,,drop=FALSE]
##	##if (existSdGeo) {
##	##	sdGEO<-matSDgeo[PTMSign$SITE,,drop=FALSE]
##	##}
##} else {
##	infoValue<-matInfoValue[PTMSign$PROTEIN,,drop=FALSE]
##	##if (existSdGeo) {
##	##	sdGEO<-matSDgeo[PTMSign$PROTEIN,,drop=FALSE]
##	##}
##}
##if (existSdGeo) {
##	write.table(sdGEO,"processed_sd_geo.txt",quote=F,sep="\t",col.name=NA)
##}
##write.table(infoValue,"processed_value_info.txt",quote=F,sep="\t",col.name=NA)
##
##
##infoValue[infoValue == "Val"]<- 2
##infoValue[infoValue == "MV"]<- 0
##infoValue[infoValue == "-Inf"]<- -1
##infoValue[infoValue == "+Inf"]<- 1
##
##write.table(infoValue,"processed_value_info_sub.txt",quote=F,sep="\t",col.name=NA)
##
##infoHM<-as.matrix(read.table("processed_value_info_sub.txt", header=T, row.names=1,check.names=1))
##
##require(pheatmap)
##png("processed_value_info.png", unit='in', width = 5, height = 15, res=600)
##pheatmap(infoHM,cluster_rows = F,cluster_cols = F, show_rownames = F, color=c("blue","black","red","white"),legend_breaks =  c(-1,0,1,2),legend_labels = c("-Infinite","Missing values","+Infinite","Values"),border_color=NA)
##dev.off()

#if (paramR$IMPUTE == "TRUE") {
#	write.table(impMat,file="processed_ratio.txt",quote=FALSE, sep="\t", col.names=NA)
#}

####>Revision history<####
# 1.1.0 [ONGOING] Major code update to remove obsolete treatments (PP 12/08/20)
# 1.0.6 Uses R_PARAM$FOCUS instead of R_PARAM$MODIF (PP 12/07/19)
# 1.0.5 sd_geo data processed only if data are available & removed/replaced 'matrix_' in files name (PP 26/03/19)
# 1.0.4 add impute data and draw quality data (SL 15/01/19)
# 1.0.3 comment AGGREGATE/GENE_NAME and MISSING_VALUES part, manage by perl script (SL 28/09/18)
# 1.0.2 Changed naming of output files (PP 26/09/18)
# 1.0.1 add drop=FALSE to keep matrix with 1 dimension (SL 24/11/17)
# 1.0.0 first version, create matrice with different option (Anova....) (SL 26/07/17)
