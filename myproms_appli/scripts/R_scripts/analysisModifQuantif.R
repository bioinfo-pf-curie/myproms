################################################################################
# analysisModifQuantif.R       1.0.3                                           #
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
#-------------------------------------------------------------------------------##Analysis Pipeline for modification and non-modification proteomics data
##No packages need

##delete previous session
rm(list=ls())
Sys.setenv( LANGUAGE = "en")

args <- commandArgs(TRUE)
matrixExport <- args[1]
pepFileExport <-args[2]

## Load R Param file from myProms
## expected format
## header :MODIF\tQUANTIF_FAM\tPROTEIN_SELECTION\tAGGREGATE\tKEEP_PROT\tKEEP_PROT_NB\tP_VALUE_CHK\tP_VALUE\tMISSING_VALUE\n
## 2ndLine:INTEGER\tRATIO\tnone(sample or group)\tTRUE(or FALSE)\tTRUE(or FALSE)\tnumber(or FALSE)\tTRUE(or FALSE)\tnumber(or FALSE)\tnumber\n

paramR<-read.table("R_parameters.txt", header=T, sep="\t")

#########Load matrix from myProms
##Matrix format
##\tsample1\tsample2\t....\tsampleN\n
##protA\tval1\tval2\t....\tvalN\n
##protB\tval1\tval2\t....\tvalN\n
#protA-1\tval1\tval2\t....\tvalN\n

matrixPTM<-read.table(matrixExport, header=T, sep="\t",row.names=1, check.names=F)

#########drop ambiguity phospho proteins (format ZFAN5_HUMAN-48~65:1/9)
##OBSOLETE MANAGE BY PERL SCRIPT
#if (paramR$EXCLUDE_AMB == "TRUE") {
#    matrixPTM<-matrixPTM[grep(":",rownames(matrixPTM),invert=T),,drop=FALSE]
#}

###############################EXPORT
#########replace infinite ratio (1000 and -1000 in this case) by NA
matrixPTM[matrixPTM == 1000]<-NA
matrixPTM[matrixPTM == -1000]<-NA

########keep only proteins with %NA in MISSING_VALUE param
#matrixPTM<-matrixPTM[rowSums(is.na(matrixPTM))<=paramR$MISSING_VALUE,,drop=FALSE]

##NONE
if (paramR$PROTEIN_SELECTION == "none") {
	PTMSign<-matrixPTM
	if (paramR$MODIF != 0){
		PTMSign$PROTEIN_MODIF<-rownames(matrixPTM)
	} else {
		PTMSign$PROTEIN<-rownames(matrixPTM)
	}
}

##SAMPLE
if (paramR$PROTEIN_SELECTION == "sample") {
    prot.sd=matrix()
    for (i in 1:length(matrixPTM[,1])) {
		prot.sd[i]=sd(matrixPTM[i,],na.rm=TRUE)
    }

    prot.sd<-as.data.frame(prot.sd)
    rownames(prot.sd)<-rownames(matrixPTM)
    colnames(prot.sd)<-"sd"
	if (paramR$MODIF != 0){
		prot.sd$PROTEIN<-gsub("-.*","", rownames(prot.sd))
		prot.sd$PROTEIN_MODIF<-rownames(prot.sd)
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
	if (paramR$MODIF != 0){
		matrixPTM<-matrixPTM[PTMSign$PROTEIN_MODIF,,drop=FALSE]
	} else {
		matrixPTM<-matrixPTM[PTMSign$PROTEIN,,drop=FALSE]
	}
}

##GROUP
if (paramR$PROTEIN_SELECTION == "group") {

    #############load group info
    ##minimum 2 samples by group
    ##the sample name in group have to be the same in matrix
    ##Group format
    ##sample\tgroup_snf
    ##sample1\tgroupA
    ##sample2\tgroupA
    ##sample3\tgroupB
    ##sample4\tgroupB
    group<-read.table("group.txt",header=T,sep="\t",row.names=1)

    ##check if there is at least 2 values by group
    nbValuesByProtBygroup<-apply(matrixPTM,1, function(x) aggregate(x, by=list(group=group$group_snf),FUN=function(y) sum(!is.na(y))>=2))

    ##get all values
    allValues<-sapply(nbValuesByProtBygroup, FUN = function(x) all(x$x) )

    ##get TRUE values
    trueValues<-as.matrix(which(allValues))

    matrixPTM<-matrixPTM[rownames(trueValues),]

    #launch ANOVA FOR GROUP
    matrixPTM_Anova<-as.data.frame(apply(matrixPTM, 1,function(x) anova(lm(x ~ group$group_snf))[1,5]))
    colnames(matrixPTM_Anova)<-"pvalue"
    if (paramR$MODIF != 0){
		matrixPTM_Anova$PROTEIN<-gsub("-.*","", rownames(matrixPTM_Anova))
		matrixPTM_Anova$PROTEIN_MODIF<-rownames(matrixPTM_Anova)
	} else {
		matrixPTM_Anova$PROTEIN<-rownames(matrixPTM_Anova)
	}
    matrixPTM_Anova<-na.omit(matrixPTM_Anova)
    nbRowAnova<-nrow(matrixPTM_Anova)

    ##P-value threshold for Anova
    if (paramR$P_VALUE_CHK == "TRUE") {
		matrixPTM_Anova<-matrixPTM_Anova[matrixPTM_Anova<paramR$P_VALUE,,drop=FALSE]
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
	if (paramR$MODIF != 0){
		matrixPTM<-matrixPTM[PTMSign$PROTEIN_MODIF,,drop=FALSE]
	} else {
		matrixPTM<-matrixPTM[PTMSign$PROTEIN,,drop=FALSE]
	}
}

#######load annotation
matrixAnnot<-read.csv("annotation.txt",sep="\t", header=T, quote = "", stringsAsFactor=FALSE)

if (paramR$MODIF != 0){##modification
	listAnnotPTM<-cbind(rownames(matrixPTM), gsub("-.*","",rownames(matrixPTM)))
	colnames(listAnnotPTM)<-c("PROTEIN_MODIF","PROTEIN")
	annotPTM<-merge(matrixAnnot, listAnnotPTM, by.x="PROTEIN", by.y="PROTEIN")
	annotPTM$PROTEIN_MODIF<-as.character(annotPTM$PROTEIN_MODIF)
} else {
	listAnnotPTM<-as.data.frame(rownames(matrixPTM))
	colnames(listAnnotPTM)<-"PROTEIN"
	annotPTM<-merge(matrixAnnot, listAnnotPTM, by.x="PROTEIN", by.y="PROTEIN")
	annotPTM$PROTEIN<-as.character(annotPTM$PROTEIN)
}

if (paramR$QUANTIF_FAM == "RATIO" && paramR$MODIF >= 0) {
	#####load matrice_pep
	matrixPepPTM<-read.table(pepFileExport,sep="\t",header=T,row.names=1,check.names=F)
	#colPep<-colnames(matrixPepPTM)
	#rowPep<-rownames(matrixPepPTM)
	#####load matrice pValue
	matrixPvaluePTM<-read.table("pvalue.txt",sep="\t",header=T,row.names=1,check.names=F)
	#colPval<-colnames(matrixPvaluePTM)
	#rowPval<-rownames(matrixPvaluePTM)

	if (paramR$MODIF != 0){
		pepPTM<-matrixPepPTM[PTMSign$PROTEIN_MODIF,,drop=FALSE]
		pValuePTM<-matrixPvaluePTM[PTMSign$PROTEIN_MODIF,,drop=FALSE]
	} else {
		pepPTM<-matrixPepPTM[PTMSign$PROTEIN,,drop=FALSE]
		pValuePTM<-matrixPvaluePTM[PTMSign$PROTEIN,,drop=FALSE]
	}

	#colnames(pepPTM)<-colPep
	#rownames(pepPTM)<-rowPep
	#colnames(pValuePTM)<-colPval
	#rownames(pValuePTM)<-rowPval
}

#if (paramR$GENE_NAME == "FALSE") {
	if (paramR$MODIF != 0) {
		annotPTM<-merge(annotPTM,PTMSign,by.x="PROTEIN_MODIF", by.y="PROTEIN_MODIF")
	} else {
		annotPTM<-merge(annotPTM,PTMSign,by.x="PROTEIN", by.y="PROTEIN")
	}
#}

##GENE_NAME, replace phospho name by gene name only if aggregate if choosen
##if (paramR$GENE_NAME == "TRUE") {
##
##	##replace empty gene by _Protein_Name
##	#annotPTM$PROTEIN_PHOSPHO<-as.character(annotPTM$PROTEIN_PHOSPHO)
##	if (paramR$MODIF != 0) {
##		annotPTM[which(annotPTM$GENE==""),"GENE"]<-annotPTM[which(annotPTM$GENE==""),"PROTEIN_MODIF"]
##		##remove duplicate GENE
##		annotPTM<-merge(annotPTM,PTMSign,by.x="PROTEIN_MODIF", by.y="PROTEIN_MODIF")
##	} else {
##		annotPTM[which(annotPTM$GENE==""),"GENE"]<-annotPTM[which(annotPTM$GENE==""),"PROTEIN"]
##		#remove duplicate GENE
##		annotPTM<-merge(annotPTM,PTMSign,by.x="PROTEIN", by.y="PROTEIN")
##	}
##
##	if (paramR$PROTEIN_SELECTION == "sample") {
##		annotPTM<-annotPTM[order(annotPTM$sd, decreasing=TRUE),]
##	}
##	if (paramR$PROTEIN_SELECTION == "group") {
##		annotPTM<-annotPTM[order(annotPTM$pvalue, decreasing=FALSE),]
##	}
##
##	annotPTM<-annotPTM[!duplicated(annotPTM$GENE),]
##
##	if (paramR$MODIF != 0) {
##		rownames(annotPTM)<-annotPTM$PROTEIN_MODIF
##		##create new matrice with no duplicated element
##		matrixPTM<-matrixPTM[annotPTM$PROTEIN_MODIF,,drop=FALSE]
##	} else {
##		rownames(annotPTM)<-annotPTM$PROTEIN
##		#create new matrice with no duplicated element
##		matrixPTM<-matrixPTM[annotPTM$PROTEIN,,drop=FALSE]
##	}
##
##	##REPLACE IN LOG MATRICE
##	rownames(matrixPTM)<-annotPTM[rownames(matrixPTM),"GENE"]
##
##	##REPLACE IN PEP MATRICE and PVALUE MATRICE
##	if (paramR$QUANTIF_FAM == "RATIO" && paramR$MODIF >= 0) {
##		if (paramR$MODIF != 0) {
##			pepPTM<-pepPTM[annotPTM$PROTEIN_MODIF,,drop=FALSE]
##			pValuePTM<-pValuePTM[annotPTM$PROTEIN_MODIF,,drop=FALSE]
##		} else {
##			pepPTM<-pepPTM[annotPTM$PROTEIN,,drop=FALSE]
##			pValuePTM<-pValuePTM[annotPTM$PROTEIN,,drop=FALSE]
##		}
##		rownames(pepPTM)<-annotPTM[rownames(pepPTM),"GENE"]
##		rownames(pValuePTM)<-annotPTM[rownames(pValuePTM),"GENE"]
##	}
##}

if (paramR$QUANTIF_FAM == "RATIO") {
	if (paramR$MODIF >= 0){
		write.table(pepPTM,paste("processed_",pepFileExport,sep=""),col.names=NA, quote=F,sep="\t")
		write.table(pValuePTM,"processed_pvalue.txt",col.names=NA, quote=F,sep="\t")
	}
	#write.table(matrixPTM,"processed_ratio.txt", col.names=NA, quote=F,sep="\t")
}
#else {
write.table(matrixPTM,paste("processed_",matrixExport,sep=""), col.names=NA, quote=F,sep="\t")
#}

if (paramR$MODIF != 0){##modification
	#annotPTM<-annotPTM[order(annotPTM$PROTEIN, decreasing=TRUE),]
	#annotPTM<-annotPTM[!duplicated(annotPTM$PROTEIN),]
	write.table(annotPTM[,c(1,2,3,4,5,6)],"processed_annotation.txt",sep="\t",quote=F, row.names = F)
} else {
	write.table(annotPTM[,c(1,2,3,4,5)],"processed_annotation.txt",sep="\t",quote=F, row.names = F)
}

if (paramR$PROTEIN_SELECTION == "sample") {
	if (paramR$MODIF != 0){##modification
		write.table(cbind(annotPTM$GENE,annotPTM$PROTEIN_MODIF,annotPTM$sd),"sd.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE","PROTEIN_MODIF","SD"))
	} else {
		write.table(cbind(annotPTM$GENE,annotPTM$PROTEIN,annotPTM$sd),"sd.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE","PROTEIN","SD"))
	}
}
if (paramR$PROTEIN_SELECTION == "group") {
	if (paramR$MODIF != 0){##modification
		write.table(cbind(annotPTM$GENE,annotPTM$PROTEIN_MODIF,annotPTM$pvalue),"anova_pvalue.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE","PROTEIN_MODIF","PVALUE"))
	} else {
		write.table(cbind(annotPTM$GENE,annotPTM$PROTEIN,annotPTM$pvalue),"anova_pvalue.txt",quote=F,sep="\t",row.names=F,col.names = c("GENE","PROTEIN","PVALUE"))
	}
}


####>Revision history<####
# 1.0.3 comment AGGREGATE/GENE_NAME and MISSING_VALUES part, manage by perl script (SL 28/09/18)
# 1.0.2 Changed naming of output files (PP 26/09/18)
# 1.0.1 add drop=FALSE to keep matrix with 1 dimension (SL 24/11/17)
# 1.0.0 first version, create matrice with different option (Anova....) (SL 26/07/17)
