################################################################################
# quantiSwath.R         2.0.14                                                 #
# Authors: Matthieu Lhotellier & Alexandre Sta (Institut Curie)                #
# Contact: myproms@curie.fr                                                    #
# Statiscal Methods for swarth protein quantification by mass spectrometry     #
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

#### Load library  ####
print("Load libraries")
library(data.table)
library(MSstats)
library(gridExtra)
library(ggplot2)
library(stringr)
library(tidyverse)
library(plyr)
library(affy)
library(reshape2)
library(FactoMineR)
library(ggdendro)
library(limma)
library(outliers)


#### SessionInfo ####
print("Save session Info")
a = toLatex(sessionInfo(), locale = FALSE)
sink(file.path("results","sessionInfo.txt"))
print(sessionInfo())
sink()

#### Load the functions ####
print("Load functions")
source(paste(filepath,"functionQuantiSwath.R",sep=""))
#### Load data  ####
print("Load Data")
# Data
areaions=fread(file.path("data","table.txt"),data.table=FALSE)

# if("excluded.txt"%in%list.files("data/")){ # DEBUG
#   areaionsExcluded=fread(file.path("data","excluded.txt"),data.table=FALSE) # DEBUG
#   areaions = areaions %>% full_join(areaionsExcluded,by=colnames(.)) # DEBUG
# } # DEBUG
#
# rowData = areaions # DEBUG

#areaions = areaions %>% as_tibble() %>% # DEBUG
#  dplyr::filter(ProteinName %in% (areaions$ProteinName %>% unique %>% .[sample(areaions$ProteinName %>% unique %>% length,100)]) ) # DEBUG

# Parameters
paramC=fread(paste(file.path("data","param_char.txt"),sep=""),header=FALSE,sep="\t",data.table=FALSE)

# PP : Make sure ProteinName is handled as character
areaions$ProteinName <- as.character(areaions$ProteinName)

# AS :  Change the ":" by "*" in the name of the peptide because MSstats can't handle that. This modification does not change the output
areaions$PeptideSequence =  gsub(":","*",areaions$PeptideSequence)
areaions$BioReplicate =  gsub("_",".",areaions$BioReplicate)  # AS :  for MissingSample()
order_names = names(areaions)

#areaions = areaions %>% filter(ProteinName=="170255") # DEBUG

#### Load Parameters ####
print("Load parameters")
source(paste(filepath,"affectParametersSwath.R",sep=""))

#### Remove peptide*fragment*protein specific to one condition across biorep&run  ####
#print("Remove peptides specific")
if(removeMissingValues){
  print("Remove peptides specific")

  # Remove the transitions specific to one condition
  cond = unique(areaions$Condition) # List of the conditions
  tmp1 = areaions[!is.na(areaions$Intensity),] # Remove all the NA from the design (the design is incomplete, it have to be full for MSstats)
  # Find the protein with missing ransitions by regrouping by dcast
  tmp1 = dcast( tmp1 ,
                ProteinName + PeptideSequence + PrecursorCharge + FragmentIon + ProductCharge + IsotopeLabelType ~ Condition ,
                value.var = "Intensity" ,
                fun.aggregate = function(x) ifelse(length(x)>0,1,0) )
  toRemove = which(apply(tmp1[,cond],MARGIN=1,FUN=sum)<2) # toRemove = the list of observations we have to remove from the data
  toRemove = tmp1[toRemove,which((!names(tmp1)%in%cond))]

  # Remove only the observations specific to one state if there is others observations common to others states in the protein
  # protSpecific = tmp1 %>% mutate(nonSpecific = State1*State2) %>% group_by(ProteinName) %>% dplyr::summarise(nonInfinity = sum(nonSpecific) ) %>% filter(nonInfinity==0)
  # toRemove = toRemove %>% filter(!ProteinName %in% protSpecific$ProteinName)

  protSpecific = tmp1 # tmp1 %>%
  protSpecific$nonSpecific = protSpecific$State1*protSpecific$State2 # mutate %>%
  protName = data.frame( ProteinName =  unique(protSpecific$ProteinName) ) # group_by %>%
  protName$nonInfinity = -1 # summarise %>%
  for(i in protName$ProteinName )  # summarise %>%
  {  # summarise %>%
    protName$nonInfinity[ which( protName$ProteinName == i ) ] = sum( protSpecific$nonSpecific[ protSpecific$ProteinName == i ] )  # summarise %>%
  }  # summarise %>%
  protName = protName[protName$nonInfinity==0,] # filter
  protSpecific = protName # filter
  toRemove = toRemove[!toRemove$ProteinName %in% protSpecific$ProteinName,] # toRemove %>%  filter


  # linesToRemove is the list of transitions which are specific to a condition
  linesToRemove = paste(toRemove$ProteinName, toRemove$PeptideSequence, toRemove$PrecursorCharge, toRemove$FragmentIon, toRemove$ProductCharge, toRemove$IsotopeLabelType )
  linesAreaions =  paste(areaions$ProteinName, areaions$PeptideSequence, areaions$PrecursorCharge, areaions$FragmentIon, areaions$ProductCharge, areaions$IsotopeLabelType )

  if(length(which(linesAreaions %in% linesToRemove))!=0){
    excluded = areaions[which(linesAreaions %in% linesToRemove),]
    areaions = areaions[which(!(linesAreaions %in% linesToRemove)),]

    # Remove protein with only NA : it come from some transitions with only NA initialy on the table.

    tmp2 = dcast(areaions,ProteinName~.,value.var = "Intensity" ,fun.aggregate =function(b) sum(!is.na(b)) )

    excluded2 = subset(areaions,ProteinName%in%tmp2$ProteinName[tmp2$.==0])
    areaions = subset(areaions,ProteinName%in%tmp2$ProteinName[tmp2$.!=0])


    # Fill the design
    #areaions$ProteinName = as.factor(areaions$ProteinName)
    #areaions$PrecursorCharge = as.factor(areaions$PrecursorCharge)
    #areaions$Run = as.factor(areaions$Run)
    #areaions$ProductCharge = as.factor(areaions$ProductCharge)

    tmp1 = dcast( areaions , ProteinName+PeptideSequence+PrecursorCharge+FragmentIon+ProductCharge+IsotopeLabelType~Condition+BioReplicate+Run, median )
    truc= melt(tmp1,id.var=1:6)
    truc$variable = as.character(truc$variable)
    truc$Condition = unlist(lapply(truc$variable,function(b)strsplit(b,split="_")[[1]][1] ) )
    truc$BioReplicate = unlist(lapply(truc$variable,function(b)strsplit(b,split="_")[[1]][2] ) )
    truc$Run = unlist(lapply(truc$variable,function(b)strsplit(b,split="_")[[1]][3] ) )
    names(truc)[which(names(truc)%in%c("variable","value"))] = c("remove","Intensity")
    truc=truc[,-which(names(truc)%in%"remove")]
    areaions = truc

    # Merge excluded & excluded2
    excluded3 = rbind(excluded,excluded2)

    # Write the data
    # write the data with the correct symbol in the sequences
    areaions.write = areaions[order(areaions$ProteinName),order_names] # order the column with the initial order & the lines by ProteinName
    areaions.write$PeptideSequence =  gsub("\\*",":",areaions.write$PeptideSequence)
    areaions.write$BioReplicate =  gsub("\\.","_",areaions.write$BioReplicate)

    excluded3.write = excluded3[,order_names] # order the column with the initial order & the lines by ProteinName
    excluded3.write$PeptideSequence =  gsub("\\*",":",excluded3.write$PeptideSequence)
    excluded3.write$BioReplicate =  gsub("\\.","_",excluded3.write$BioReplicate)

    write.table(areaions.write,file.path("data","table.txt"),quote=FALSE,row.names=FALSE,sep="\t")
    write.table(excluded3.write,file.path("data","excluded.txt"),quote=FALSE,row.names=FALSE,sep="\t")
  }
}
#
# nprot = 250 # DEBUG
# n = length(unique(areaions$ProteinName)) # DEBUG
# areaions = areaions %>% filter(ProteinName%in% c(unique(areaions$ProteinName)[sample(1:n,nprot)],"3583191") ) # DEBUG

#### Check Data type ####

#### Description Raw Data  ####
print("Description of the data")
data.state="brut"
normalisation = "Before"
visualisation = TRUE
missingsample = TRUE

# Script
dataVis = areaions
source(paste(filepath,"visualisation.R",sep=""))


#### Normalisation & description of the data  ####

print("Normalization")
QuantData<-dataProcess( areaions,
                        normalization = normalization,
                        nameStandards = protNorm,
                        clusters = as.numeric(clusters),
                        featureSubset = featureSubset,
                        n_top_feature = n_top_feature,
                        summaryMethod = summaryMethod
                        )

#### Visualisation
# parameters
data.state="After"
normalisation = "After"
visualisation = TRUE
missingsample = FALSE
# script
dataVis = QuantData$ProcessedData %>% select(PROTEIN,PEPTIDE,TRANSITION,FEATURE,LABEL,GROUP_ORIGINAL,SUBJECT_ORIGINAL,RUN,ABUNDANCE) %>%
  separate(TRANSITION,c("FragmentIon","ProductCharge"),sep="_") %>%
  separate(PEPTIDE,c("PeptideSequence","PrecursorCharge"),sep="_") %>%
  dplyr::rename(ProteinName=PROTEIN,IsotopeLabelType=LABEL,Intensity=ABUNDANCE,Condition=GROUP_ORIGINAL,BioReplicate=SUBJECT_ORIGINAL,Run=RUN) %>%
  mutate(Intensity=2^Intensity)
source(paste(filepath,"visualisation.R",sep=""))

#### Outlier ####
print("Outlier")
outlier_detect = FALSE # Dont Remove the outliers !!!
if(outlier_detect){
  QuantData<-HandlingCVout(QuantData)
}

#### Visualisation
# parameters
#data.state="Afterout"
#normalisation = "After_CVout"
#visualisation = FALSE
#missingsample = FALSE
# script
#source(paste(filepath,"visualisation.R",sep=""))



#### Remove prot one cond ####
#QuantData=removeonecond(QuantData,level="prot") # Remove lines of prot from the data where prot is unique to only one condition

#Export Pep used for quantification
#exportpep(QuantData)


#### Build the contrast matrix  ####
print("Build contrast matrix")

if( length(matrix.contrast) == 0 ){
  groups = levels(QuantData$ProcessedData$GROUP_ORIGINAL)
  grps = groups[sapply(X = groups,.grepNum) %>% as.numeric() %>% order(decreasing=TRUE)] %>% combn(2) %>% t %>% as_tibble %>%
    unite(names,c("V1","V2"),sep="-",remove=FALSE)

  contrast.matrix = matrix(0,ncol=length(unique(groups)),nrow=dim(grps)[1])
  rownames(contrast.matrix) = grps$names
  colnames(contrast.matrix) = groups
  for(i in 1:dim(contrast.matrix)[1] ){
    left = which(colnames(contrast.matrix)==unlist(strsplit(rownames(contrast.matrix)[i],split="-"))[1])
    right = which(colnames(contrast.matrix)==unlist(strsplit(rownames(contrast.matrix)[i],split="-"))[2])
    contrast.matrix[i,left] = 1
    contrast.matrix[i,right] = -1
  }
}else{
  contrastRowName = matrix.contrast %>% strsplit(split=";") %>% unlist() %>% gsub("/","-",.)
  contrastColName = contrastRowName %>% paste(collapse="-") %>% strsplit(split="-") %>% unlist %>% unique
  contrast.matrix = matrix( 0,nrow = length(contrastRowName) , ncol = length(contrastColName) )
  rownames(contrast.matrix) = contrastRowName
  colnames(contrast.matrix) = contrastColName
  for(i in contrastRowName){
    LG = strsplit(i,split="-") %>% unlist
    contrast.matrix[i,LG[1]] = 1
    contrast.matrix[i,LG[2]] = -1
  }
}

# Force the order of the columns to be alphabetically ordered (specified in groupComparison)
tmpRowNames = rownames(contrast.matrix)
tmpColNames = colnames(contrast.matrix)
tmpDim = dim(contrast.matrix)
contrast.matrix = contrast.matrix[,order(contrast.matrix %>% colnames())]
if(class(contrast.matrix)=="numeric"){
  contrast.matrix = matrix(contrast.matrix,tmpDim)
  rownames(contrast.matrix) = tmpRowNames
  colnames(contrast.matrix) = tmpColNames
}else{
  rownames(contrast.matrix) = tmpRowNames
}


#### Differential Analysis  ####
print("AnalysisDiff") # AS
res = groupComparison( contrast.matrix = contrast.matrix , data = QuantData )

#### Change the NA in Inf ####
# AS character for speeding the code
res$ComparisonResult$Protein = as.character(res$ComparisonResult$Protein)
res$QuantData$PROTEIN = as.character(res$QuantData$PROTEIN)

# Tronc the dataset to launch the code only on a subset
infOrNa = res$ComparisonResult %>% filter(is.na(log2FC))
QD_tronc = QuantData$ProcessedData %>% filter(PROTEIN %in% infOrNa$Protein)

if(dim(infOrNa)[1]!=0){
  # Replace the NA by the appropriate information (Inf, -Inf or NA)
  fun<-function(b){
    lab = unlist(strsplit(as.character(b[,"Label"]),split="-"))
    tmp = subset(QD_tronc,PROTEIN==b[,"Protein"]&GROUP_ORIGINAL%in%lab&!is.na(ABUNDANCE))
    pos = dim(subset(tmp,GROUP_ORIGINAL==lab[1]))[1]
    neg = dim(subset(tmp,GROUP_ORIGINAL==lab[2]))[1]

    if(pos*neg>0){ # When there peptides or transitions do not match between the 2 conditions, just return NA
     return(b)
    }
    b[,"log2FC"]=ifelse((pos-neg)==0,NA,sign(pos-neg)*Inf)
    return(b)
  }

  infOrNa=adply(.data=infOrNa,.margin=1,.fun=fun,.progress="text")


  # Fill the original matrix ComparisonResluts with the new correct values of log2FC

  for(i in 1:dim(infOrNa)[1])
  {
    ind = which(res$ComparisonResult$Protein==infOrNa$Protein[i] & res$ComparisonResult$Label==infOrNa$Label[i])
    res$ComparisonResult[ind,"log2FC"] = infOrNa[i,"log2FC"]
  }
}

res$ComparisonResult$Protein = as.factor(res$ComparisonResult$Protein)

#### Add the Inf missing from res ####
# List of all the protein Inf
listOfInf = areaions %>% group_by(ProteinName) %>% nest %>% mutate( Label=map(1,function(x)rownames(contrast.matrix)) ) %>% select(-data) %>% unnest(Label) %>%
  full_join((areaions %>% group_by(ProteinName) %>% nest),by="ProteinName") %>% dplyr::rename(Protein=ProteinName) %>%
  separate(Label,c("left","right"),remove=FALSE) %>% mutate( conditions = map(data,function(x) unique(x %>% filter(!is.na(Intensity)) %>% .$Condition)) ) %>%
  select(-data) %>% mutate(conditions=map(conditions,function(x)paste(x,collapse="."))) %>%
  mutate( conditions = unlist(conditions)) %>% unite(leftCond,c("left","conditions") , sep="/" , remove=FALSE ) %>%
  unite( rightCond , c("right","conditions") , sep="/" , remove=FALSE ) %>% select( -right , -left , -conditions ) %>%
  mutate( leftCond = map( leftCond,function(x) strsplit(x,split="/")[[1]][1]%in%unlist(strsplit(strsplit(x,split="/")[[1]][2],split = "\\.")) ) ) %>%
  mutate( rightCond = map( rightCond,function(x) strsplit(x,split="/")[[1]][1]%in%unlist(strsplit(strsplit(x,split="/")[[1]][2],split = "\\.")) ) ) %>%
  mutate( leftCond = unlist(leftCond) , rightCond = unlist(rightCond) ) %>% filter(!leftCond|!rightCond) %>%
  mutate( leftCond = ifelse(leftCond,1,0) , rightCond = ifelse(rightCond,1,0) ) %>% mutate(log2FC=log(leftCond/rightCond)) %>%
  select(-rightCond,-leftCond) %>% mutate( issue = "oneConditionMissing_Added" )

# Select only protein specific to listOfInf and add them to res$ComparisonResult
#res$ComparisonResult = listOfInf %>% anti_join( (res$ComparisonResult %>% mutate(Protein=as.numeric(as.character(Protein))) %>% select(Protein,Label)) , by=c("Protein","Label") ) %>%
#  full_join(res$ComparisonResult %>% mutate(Protein=as.numeric(as.character(Protein)),Label=as.character(Label),issue=as.character(issue)) , by=c("Protein","Label","log2FC","issue") )
res$ComparisonResult = listOfInf %>% anti_join( (res$ComparisonResult %>% mutate(Protein=as.character(Protein)) %>% select(Protein,Label)) , by=c("Protein","Label") ) %>%
  full_join(res$ComparisonResult %>% mutate(Protein=as.character(Protein),Label=as.character(Label),issue=as.character(issue)) , by=c("Protein","Label","log2FC","issue") )

#### Export the final results ####
# Replace the NaN with NA
res$ComparisonResult$log2FC[is.nan(res$ComparisonResult$log2FC)] = NA
# Filter the NA log2FC in the results, those proteins will not be reported.
res$ComparisonResult = res$ComparisonResult %>% filter(!is.na(log2FC))

AD = res$ComparisonResult %>% select(Protein,Label,log2FC) %>% dplyr::rename(log2=log2FC,log2FC=Label) %>%  spread(log2FC,log2,sep = "_") %>%
  full_join( (res$ComparisonResult %>% select(Protein,Label,pvalue) %>% dplyr::rename(pval=Label) %>%  spread(pval,pvalue,sep = "_")) , by="Protein" ) %>%
  full_join( (res$ComparisonResult %>% select(Protein,Label,SE) %>% dplyr::rename(SE_=SE,SD=Label) %>%  spread(SD,SE_,sep = "_")) , by="Protein")

if(sum(! (is.nan(res$ComparisonResult$pvalue)|is.na(res$ComparisonResult$pvalue))) ){
  groupComparisonPlots(data=res$ComparisonResult,type="VolcanoPlot",logBase.pvalue=10,address="results/graph/",ProteinName = FALSE)
  p = res$ComparisonResult %>% ggplot(aes(adj.pvalue)) + geom_histogram(binwidth=0.01) + geom_density() +  facet_wrap(~Label)
  ggsave("results/graph/pvalHist.jpeg",width=10,height=10)
}

write.table(AD,file.path("results","ResultsDAProt.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(res$ComparisonResult,file.path("results","result1.txt"),quote=FALSE,row.names=FALSE,sep="\t")

print("END")

####>Revision history<####
# 2.0.14 Force "ProteinName" column type to character to prevent default attribution as numeric (PP 18/07/18)
# 2.0.13 little bug removed from the last modification (AS 01/03/18)
# 2.0.12 correction of the filte peptide*fragment*protein specific for multiple states (AS 01/03/18)
# 2.0.11 remove debug lines (AS 09/02/18)
# 2.0.10 remove the non quantified proteins (AS 09/02/18)
# 2.0.9 when there is no peptides matching between conditions, put NA in the data (AS 09/02/18)
# 2.0.8 works with both two samples or more than two conditions (AS 08/02/18)
# 2.0.7 contrast matrix is now ordered alphabetically (AS 07/02/18)
# 2.0.6 the version 2.0.5 without the bug (AS 26/01/18)
# 2.0.5 correct the construction of the normalization matrix and speed up the visualisation step (AS 24/01/18)
# 2.0.4 add options in dataProcess and debug the cluster option (AS 19/01/18)
# 2.0.3 add Inf in protein specific missing from the data (AS 16/01/18)
# 2.0.2 add the cluster option for MSstats (AS 15/01/18)
# 2.0.1 when there is no bioReplicate, do not repport a bug (AS 22/11/17)
# 2.0.0 the code is now compatible with MSstats >3 (AS 21/11/17)
# 1.2.7 remove data.table:: before dcast (AS 19/09/17)
# 1.2.6 addapt the code when the library tidyr and dplyr is not allowed (AS 19/09/17)
# 1.2.5 change the correction method of the p.value to fdr (AS 06/09/17)
# 1.2.4 Correction on the fragment specific filter (AS 23/03/17)
# 1.2.3 Remove the file CVout.txt (AS 15/12/16)
# 1.2.2 Reorder the lines of the tables table.txt & excluded.txt (AS 08/12/16)
# 1.2.1 Reorder the column of the tables table.txt & excluded.txt (AS 07/12/16)
# 1.2.0 Correct the symbols for the modifications in the peptide sequence for output table (AS 07/12/16)
# 1.1.9 Commentar and minor modification of the code for the step : remove outliers (AS 06/12/16)
# 1.1.8 Merge excluded & excluded 2 to have onl one file excluded (AS 05/12/16)
# 1.1.7 Fix the bug which erase table.txt, now it correctly removes the specific fragments (AS 01/12/16)
# 1.1.6 Fix the bug log2FC, constraint the name of this column to be log2FC and not logFC (AS 29/11/16)
# 1.1.5 exclusion of the missing transitions reported in the datafile : excluded.txt (AS 18/11/16)
# 1.1.4 globalStandards normalization is aviable (AS 12/08/16)
# 1.1.3 normalize correctly & print time in .rout (AS 04/08/16)
# 1.1.2 in ResultDAProt.txt don't remove protein specific to one condition but add Inf, -Inf or NAs (AS 03/08/16)
# 1.1.1 add default script for matrix.contrasts (AS 02/08/16)
# 1.1.0 change _ in bioRep, change the name convention of the script, add scripts & read paramC (AS 01/08/16)
# 1.0.4 modification of the names of the output files & remove : from data (AS 27/06/16)
# 1.0.3 modification of filepath and remove 00- & 01- (AS 18/06/16)
# 1.0.2 correction filepath for outputs & add an output of sessionInfo (AS 06/06/16)
# 1.0.1 correction of the normalisation parameter : constant  (AS 04/11/15)
