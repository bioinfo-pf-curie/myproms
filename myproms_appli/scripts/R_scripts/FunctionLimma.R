################################################################################
# FunctionLimma.R         5.2.4                                              #
# Authors: Matthieu Lhotellier & Alexandre Sta & Isabel Brito (Institut Curie) #
# Contact: myproms@curie.fr                                                    #
# Function of AnalysisDiffLimma.R                                              #
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

#####
# Load Packages
# require(affy) # Function normalize.loess, commented because it is too long to load. I build my own loess normalization function
library(plyr) # Function ddply, think to replace the function ddply by a new one from tidyverse to remove this package
library(ggplot2) # To build beautifull plots
library(tidyr) # To reshape the data with powerfull & readable function
library(data.table) # To read the data quicklier
library(purrr) # To have the map function
library(dplyr) # Function full_join
library(preprocessCore) # Function normalizeQuantile
library(broom) # To make tidy a linear model from lm (glance)
library(stringr) # To manipulate character strings
library(igraph) # To rebuild the correct names of the replicates
library(readr) # for the function ???
library(bit) # for fread  
library(bit64) # for fread 

#################################################################################
## Name : .myFilterOutlier
## Arguments   : object
## Description : 
#################################################################################
.myFilterOutlier = function(object,params){
  print("Execution of .myFilterOutlier")
  X = list(data = object , outlier = object %>% filter(proteinId==-1) )
  
  # Remove the biological outliers
  #X = .cvOutliers(X,rep="bio")
  
  # Remove the technical outliers
  #X = .cvOutliers(X,rep="tech")
  # Remove the peptide outliers
  
  switch(params$design,
         PEP_RATIO = {
           X = .peptideOutliers(X)
         },
         PEP_INTENSITY = {
           X = .mixBioPeptideOutliers(X,params)
         }
  )
  
  # Remove the missing sample
  # if(parameters$design!="PEP_RATIO")
  # {
  #   X$data = .removePeptideMissing(X$data)
  # }
  
  print(".myFilterOutlier done")
  return(X)
}

#################################################################################
## Name : .reshapeData
## Arguments   : path
##               print
## Description : 
#################################################################################
.reshapeData <- function(object)
{
  # Select only the usefull columns :
  object = object %>%
    select("Protein_ID","Peptide","Sample","Replicate","Technical_Replicate","Experiment","Quantif_Set","Peptide_IDs","Value")
  
  # Add the column Experiment if it does not exist
  if( !( "Experiment" %in% names( object ) ) | !sum( !is.na( object$Experiment ) ) ){ object$Experiment = "A" }
  # Rename the columns with the convention word1Word2Word3...
  orderData = data.frame(
    key = names( object )
  )
  orderRename = data.frame(
    key = c("Peptide","Sample","Protein_ID","Replicate","Technical_Replicate","Experiment","Quantif_Set","Peptide_IDs","Value"),
    namesNew = c("peptide","sample","proteinId","replicate","repTech","experiment","quantifSet","peptideId","value")
  )
  fullJoin = full_join( orderData , orderRename , by = "key" )
  if( sum( is.na( fullJoin ) ) ){ stop("wrong names or number of columns") }
  names( object ) = fullJoin$namesNew
  
  # Remove the characters after the last "_" in the column peptide
  object = object %>% separate( peptide , c( "peptide" , "first" , "second" ) , sep = "_" ) %>%
    unite( peptide , peptide , first , sep = "_" ) %>%
    dplyr::select( -second )
  # Select only the good columns (in particular remove treatment)
  object = object %>% select( proteinId,
                              peptide,
                              sample,
                              replicate,
                              repTech,
                              experiment,
                              quantifSet,
                              peptideId,
                              value
  )
  return(object)
}

#################################################################################
## Name : .saveSessionInfo
## Arguments   : path
##               print
## Description : 
#################################################################################
.saveSessionInfo <- function( path = "" , print = TRUE )
{
  sink( file.path( path , "sessionInfo.txt" ) )
  print( sessionInfo( ) )
  sink( )
  
  if( print )
  {
    print( sessionInfo( ) )
  }
}

#################################################################################
## Name : .recoverReplicate
## Arguments   :
##               "dataTmp" type data.frame
## Description : Rebuild the number of the run from the columns quantifSet, sample
##               replicate and repTech
##
#################################################################################
.recoverReplicate = function(dataTmp){
  # Build the matrix of relations
  tmp  = dataTmp %>% select(sample,replicate,quantifSet) %>% unite(rep,sample,replicate,remove=FALSE,sep="") %>%
    select(-replicate) %>%
    spread(key = sample,value=rep) %>% select(-quantifSet) %>% unique
  
  # Reshape the matrix of relations to have a graph
  aa = tmp %>% names %>% combn(2) %>% t
  run = data_frame( name = unique( unlist(tmp)[!is.na(unlist(tmp))] ) )
  relations = data_frame(from=character(0),to=character(0))
  for(i in 1:dim(aa)[1])
  {
    relTmp = tmp[,aa[i,]] %>% na.omit  ; names(relTmp) = c("from","to")
    relations = relations %>%
      rbind(relTmp)
  }
  g <- graph_from_data_frame(relations, directed=FALSE, vertices=run)
  # Find the connex part of the graph
  clust =  clusters(g)$membership
  dataTmp = dataTmp %>% unite(rep,sample,replicate,sep="",remove=FALSE) %>% select(-replicate) %>%
    mutate( replicate = clust[rep] ) %>% select(-rep)
  return(dataTmp)
}


#################################################################################
## Name : .myUnSplit
## Arguments   :
##               "x"
## Description :
##
#################################################################################
.myUnsplit <- function(x)
{
  v = seq(1,length(x),2500)
  left = v[1:(length(v)-1)]
  right = v[2:length(v)]+1
  right[length(v)-1] = tail(right,1)-1
  
  res_tmp = list(NULL)
  for(i in length(left))
  {
    res_tmp[[i]] = do.call("rbind", x[left[i]:right[i]]) # store the result in one data frame
  }
  return(do.call("rbind", res_tmp))
}
#################################################################################
## Name : .myRename
## Arguments:
##           ".data" tibble or data.frame, the data
##           "oldNames" character , names of the columns of data you want to change
##           "newNames" character, new names you want (oldNames & newNames have to
##                      be the same lenght)
## Description :
#################################################################################
.myRename_ <- function(.data,oldNames,newNames){
  out <- dplyr::rename_(.data, .dots=setNames(oldNames, newNames))
  return(out)
}

#################################################################################
## Name : .missingSample
## Arguments:
##           "object" data
##           "step" step of workflow (ex: "Before_Normalization")
## Description :
#################################################################################
.missingSample <- function(object,step){
  print("Execution of the function .missingSample")
  # Build the table where each lines are a unique : protein+peptide+sample+replicate+repTech+experiment,
  # all the redondant values of log2-intensity are merged in value with the function mean
  tab=object
  tab = tab %>% group_by(proteinId,peptide,sample,replicate,repTech,experiment) %>%
    dplyr::summarise(value = mean(value,na.rm=TRUE)) %>% ungroup %>% setDT()
  #tab = tab %>%
  #  aggregate(value~proteinId+peptide+sample+replicate+repTech+experiment, data=., FUN=mean) %>% # merge the multiple obsevation per (formula in aggregate)
  #  setDT()  # Change the data into a table by reference
  # Reshape the data where each lines are proteinId+peptide, the columns are the experiment/sample/replicate/repTech
  res = tab %>%
    dcast.data.table(proteinId+peptide~experiment+sample+replicate+repTech,value.var="value",sep="/") %>%
    setDF()  # Change the data into a data.frame by reference
  # Write the percentage of NA in percentageNAsample.txt
  x = res %>% select(-proteinId,-peptide) %>% apply(2,function(y)return(sum(is.na(y))))
  write.table( x*100/dim(res)[1],
               file.path("results",paste(step,"percentageNAsample.txt",sep="")),
               quote=FALSE,
               row.names=TRUE,
               col.names=FALSE,
               sep="\t"
  )
  # Build the table for the image
  tmp = data.frame(sample = names(x*100/dim(res)[1]) ,value = x*100/dim(res)[1])
  tmp = separate(tmp,sample,c("experiment","state","replicate","techRep"),sep="/",remove=FALSE)
  # Plot and save
  ggplot(tmp,aes(x = sample,y=value,fill=state)) + geom_bar(stat = "identity") + ylab("Percentage of NA")
  ggsave(file=file.path("results","graph",paste(step,"nbrowmissing.jpeg",sep="")),width=11,height=9)
  graphics.off()
}

#################################################################################
## Name : .calculLog2Ratio
## Arguments:
##           "object" data
##           "parameters" the list of the parameters
## Description : Compute all the possible ratios & remove the column treatment
#################################################################################
.calculLog2Ratio<-function(object,parameters){
  print("Execution of .calculLog2Ratio")
  design = parameters$design
  
  if(design!="PEP_RATIO"){
    # Aggregation step : aggregate the multiple values of the a peptide (due to the multiple integration of the spectrum)
    #                   by adding all the intensity values observed (addition on the initial intensity not on the log2)
    object = object %>%  # reverse the log2 transformation
      group_by(proteinId,peptide,experiment,replicate,repTech,sample) %>%
      dplyr::summarise(value=sum(as.numeric(value)), # HERE TO CHANGE THE FUNCTION OF AGGREGATION, as.numeric is to convert the data to numeric because the integer in R are coded in 32 bits
                       quantifSet = paste(quantifSet,collapse="+"),
                       peptideId = paste(peptideId,collapse="+")
      ) %>% mutate(value = log2(value)) %>% ungroup # aggregation and log2 transformation
  }
  
  switch(design,
         PEP_RATIO = {
           print("Design PEP_RATIO : build the ratio of the observations")
           # # Long version of the computation of the ratios (grouping by Qset and compute the pairs)
           # 
           # fct <-function(x){
           #   ref = x %>% filter(sample=="State1")
           #   x = x %>% filter(sample!="State1") %>% mutate(M=value-ref$value,
           #                                             A=(value+ref$value)/2,
           #                                             sample=paste(sample,ref$sample,sep="."),
           #                                             peptideId=paste(peptideId,ref$peptideId,sep=".")
           #   )
           #   return(x)
           # }
           # object = object %>% group_by(quantifSet) %>% nest %>%
           #   mutate(n = map(data,function(x)dim(x)[[1]])) %>% unnest(n) %>% filter(n!=1) %>% select(-n) %>%
           #   mutate(others=map(data,fct)) %>% select(-data) %>% unnest %>%
           #   select(-value)
           
           
           object = object %>% mutate(value = log2(value))
           # Short version of computing the ratios (build a ref data table and the others observations and merge by Qset)
           others = object %>% filter(sample!="State1")
           ref = object %>% filter(sample=="State1")
           object = left_join(others,ref,by=c("quantifSet","peptide","proteinId","experiment")) %>%
             select(-replicate.y,-repTech.y) %>% unite(sample,sample.x,sample.y,sep=".") %>%
             mutate( M=value.x-value.y ,
                     A=(value.x+value.y)/2 ,
                     replicate=replicate.x,
                     repTech = repTech.x,
                     peptideId = paste(peptideId.x,peptideId.y,sep=".")
             ) %>%
             select(-value.x,-value.y,-replicate.x,-repTech.x,-peptideId.x,-peptideId.y) %>%
             filter( !is.na(M))
           
           # Aggregation step : aggregate the multiple values of the a peptide (due to the multiple integration of the spectrum)
           object = object %>%
             group_by(proteinId,peptide,experiment,replicate,repTech,sample) %>%
             dplyr::summarise(M = median(M), # HERE TO CHANGE THE FUNCTION OF AGGREGATION
                              A = median(A),
                              quantifSet = paste(quantifSet,collapse="+"),
                              peptideId = paste(peptideId,collapse="+")
             ) %>% ungroup 
           
           print(".calculLog2Ratio done")
           return(object)
         },
         PEP_INTENSITY = {# In PEP_INTENSITY with multiple reference an observation is an intensity
           print("Design PEP_INTENSITY : don't compute a ratio")
           object = object %>% ungroup %>%
             mutate(A = 0,M=value) %>% select(-value)
           print(".calculLog2Ratio done")
           return(object)
         }
  )
}


#################################################################################
## Name : .removePeptideMissing
## Arguments:
##           "object" data
## Description : Remove the peptides missing in allmost one conditon in all the
##               the others conditions
#################################################################################
.removePeptideMissing <- function(object){
  
  samp = object %>% select(experiment,sample) %>% unique %>% group_by(experiment) %>% dplyr::summarise(nThSample=n())
  
  object = object %>% group_by(proteinId,peptide,experiment) %>% nest %>% 
    mutate( nObsSample = map(data,function(x) length(unique(x$sample)) ) ) %>% 
    unnest( nObsSample ) %>% full_join(samp,by="experiment") %>% 
    mutate(remove = nObsSample-nThSample) %>% filter(remove==0) %>% select(-remove,nObsSample,nThSample) %>% 
    unnest(data)
  
  return(object)
}

#################################################################################
## Name : .normalizeData
## Arguments:
##           "object" data.frame
##           "params" data.frame
##           "normProt" character
## Description : Normalise the runs one to each others
#################################################################################
.normalizeData <- function(object,params,normProt){
  print("Execution of .normalizeData")
  # Read the parameters of the normalisation  
  normalizationMethod = unlist(strsplit(params$normalization.method,split="\\."))
  if(length(normalizationMethod)==1){normalizationMethod=c(normalizationMethod,"none")}
  # If normProtein is NULL use all the proteins to normalise
  if(is.null(normProt)){
    normProtein=data.frame(proteinId=unique(object$proteinId))
  }
  
  # Tag the specific proteins who will be used to normalize the data (possibly all, at least one)
  object$normProtein = 0
  object$normProtein[ which( object$proteinId %in% normProtein$proteinId ) ] = 1
  
  # Select the functions to normalize
  if(normalizationMethod[1] != "quantile"){
    switch( normalizationMethod[1],
            none = {
              E = function(x) return(0)
            },
            median = {
              E = median
            },
            mean = {
              E = mean
            }
    )
    switch(normalizationMethod[2],
           none = {
             V = function(x) return(1)
           },
           scale = {
             V = mad
           }
    )
    fctNorm = .myNormWB
    
  }else{
    fctNorm = .myQuantileNorm
  }
  
  # Start Normalisation
  print("Normalisation processing ...")
  resu = fctNorm(object,E,V,design=params$design)
  
  # # Start DEBUG section
  # data$M[data$experiment=="A"] = data$M[data$experiment=="A"] + 
  #   rnorm(sum(data$experiment=="A"),2,0.3)
  # data$M[data$experiment=="A"&data$replicate=="rep2"] = data$M[data$experiment=="A"&data$replicate=="rep2"] + 
  #   rnorm(sum(data$experiment=="A"&data$replicate=="rep2"),1,1)
  # data$M[data$experiment=="A"&data$replicate=="rep3"] = data$M[data$experiment=="A"&data$replicate=="rep3"] + 
  #   rnorm(sum(data$experiment=="A"&data$replicate=="rep3"),-1,0.3)
  # data$M[data$experiment=="B"&data$replicate=="rep2"] = data$M[data$experiment=="B"&data$replicate=="rep2"] + 
  #   rnorm(sum(data$experiment=="B"&data$replicate=="rep2"),1.5,2)
  # data$M[data$experiment=="B"&data$replicate=="rep2"] = data$M[data$experiment=="B"&data$replicate=="rep2"] + 
  #   rnorm(sum(data$experiment=="B"&data$replicate=="rep2"),4,3)
  # 
  # library(cowplot) # DEBUG
  # p1 = data %>% unite(x,replicate,sample,remove = FALSE) %>%  ggplot(aes(x,M,color=sample)) + geom_boxplot() +
  #   facet_wrap(~experiment)+ ggtitle("Row Data") + xlab("sample") + ylab("log2(intensity)")  # DEBUG
  # p2 = tmpAdd$data %>% unite(x,replicate,sample,remove = FALSE) %>%  ggplot(aes(x,M,color=sample)) + geom_boxplot() +
  #   facet_wrap(~experiment) + ggtitle("Additive bias correction") + xlab("sample") + ylab("log2(intensity)") # DEBUG
  # p3 = tmpMul$data %>% unite(x,replicate,sample,remove = FALSE) %>%  ggplot(aes(x,M,color=sample)) + geom_boxplot() +
  #   facet_wrap(~experiment) + ggtitle("Additive & multiplicative bias correction") + xlab("sample") + ylab("log2(intensity)")# DEBUG
  # plot_grid(p1,p2,p3,ncol = 1) # DEBUG
  # historyData$dataRatio %>% full_join( # DEBUG
  #   (bias %>% filter(!run%in%LETTERS) %>% separate(run,c("experiment","sample","replicate"))), # DEBUG
  #   by=c("experiment","sample","replicate") # DEBUG
  #   ) %>% mutate(Mprime=coefMul*M+coefAdd) %>% group_by(sample,replicate) %>% # DEBUG
  #   dplyr::summarise(m=median(Mprime),sd=mad(Mprime)) # DEBUG
  # # End DEGUG section
  print(".normalizeData done")
  return(resu)
}

#################################################################################
## Name : .myNormWB 
## Arguments:
## Description :
#################################################################################
.myNormWB <-function(object,E,V,design){
  # Normalize within experiment 
  X = .normWithinExp(object,E,V,design)
  object = X$data
  bias = X$bias
  # Normalize between experiment
  X = .normBetweenExp(object,E,V,design)
  object = X$data
  bias = bias %>% unite(run,c("experiment","sample","replicate","repTech"),sep="_") %>% full_join( (X$bias %>% dplyr::rename(run=experiment) ) , by=names(.) )
  
  return(list(data=object,bias=bias))
}

#################################################################################
## Name : .normBetweenExp 
## Arguments:
## Description :
#################################################################################
.normBetweenExp <- function(object,E,V,design){
  # Build m_i and \sigma_i
  biasExp = object %>% filter(normProtein==1) %>% group_by(experiment) %>% 
    dplyr::summarise(coefMul = V(M),coefAdd=E(M))
  biasRef = biasExp %>% group_by() %>% dplyr::summarise(coefMulRef = exp(mean(log(coefMul))) , coefAddRef = mean(coefAdd) ) 
  if( (design=="PEP_RATIO") & (E(c(1,2,3))!=0) ){ 
    biasRef$coefAddRef = 0
  } 
  bias = biasExp %>% mutate(by=1) %>%  full_join( (biasRef %>% mutate(by=1)), by="by" ) %>%
    mutate(coefAdd=coefAddRef-(coefMulRef/coefMul)*coefAdd,coefMul=coefMulRef/coefMul) %>%
    select(-coefMulRef,-coefAddRef,-by)
  bias.output = biasExp %>% mutate(by=1) %>%  full_join( (biasRef %>% mutate(by=1)), by="by" ) %>%
    mutate(coefAdd=coefAddRef-coefAdd,coefMul=coefMulRef/coefMul) %>%
    select(-coefMulRef,-coefAddRef,-by)
  # Add the real coefficients of correction of the normalization
  bias.output = bias.output %>% mutate(by=1) %>% full_join( (bias %>% mutate(by=1) %>% select(-coefMul) %>% dplyr::rename(coefAddReal="coefAdd")) , by = c("by","experiment")  ) %>% 
    select(-by)
  # Correct and build the output
  object = object %>% full_join(bias,by=c("experiment")) %>% mutate(M=coefMul*M+coefAdd) %>% select(-coefMul,-coefAdd)
  X = list(data=object,bias=bias.output)
  return(X)
}

#################################################################################
## Name : .normWithinExp 
## Arguments:
## Description :
#################################################################################
.normWithinExp <- function(object,E,V,design){
  # Build m_i and \sigma_i
  biasRun = object %>% filter(normProtein==1) %>% group_by(sample,replicate,experiment,repTech) %>% 
    dplyr::summarise(coefMul = V(M),coefAdd=E(M))
  biasExperiment = biasRun %>% group_by(experiment) %>% dplyr::summarise(coefMulExp = exp(mean(log(coefMul))) , coefAddExp = mean(coefAdd) )
  # If PEP_INTENSITY center the mean to the mean of the mean, if PEP_RATIO m_i = 0
  if( (design=="PEP_RATIO") & (E(c(1,2,3))!=0) ){ 
    biasExperiment$coefAddExp = 0
  } 
  bias.output = biasRun %>% full_join(biasExperiment,by="experiment") %>%
    mutate(coefAdd=coefAddExp-coefAdd,coefMul=coefMulExp/coefMul) %>%
    select(-coefMulExp,-coefAddExp) 
  bias = biasRun %>% full_join(biasExperiment,by="experiment") %>%
    mutate(coefAdd=coefAddExp-(coefMulExp/coefMul)*coefAdd,coefMul=coefMulExp/coefMul) %>%
    select(-coefMulExp,-coefAddExp)
  # Add the real coefficients of correction of the normalization
  bias.output = bias.output %>% full_join( (bias %>% select(-coefMul) %>% dplyr::rename(coefAddReal="coefAdd")) )
  # Build the output
  object = object %>% full_join(bias,by=c("sample","replicate","experiment","repTech")) %>% mutate(M=coefMul*M+coefAdd) %>% select(-coefMul,-coefAdd)
  X = list(data=object,bias=bias.output)
  return(X)
}



#################################################################################
## Name : ..myQuantileNorm
## Arguments:
## Description :
#################################################################################
.myQuantileNorm=function(object,...){
  objectSave = object
  # Build a matrix with protein//peptide in line and the sample in column for the normalization
  object = object %>%
    select(proteinId,peptide,M,experiment,replicate,repTech,sample) %>%
    unite(run,experiment,replicate,repTech,sample,sep="//") %>%
    spread(run,M)
  # Normalize
  tmp = object %>% select(-proteinId,-peptide) %>% as.matrix %>%
    normalize.quantiles(copy=FALSE)  
  object = cbind( object[,c("proteinId","peptide")] ,tmp) %>% as_tibble %>% 
    gather_(key_col = "sample" , value_col = "M",colnames(tmp)) %>% 
    separate( sample , c("experiment","replicate","repTech","sample") , sep = "//" ) %>%
    dplyr::right_join( ( objectSave %>% select(-M) ) , by = c("proteinId","peptide","experiment","replicate","repTech","sample") )
  resu = list(data=object,bias=NULL)
  return(resu)
}

#################################################################################
## Name : .globModel
## Arguments:
##           "object" data.frame
##           "parameters" data.frame
## Description :
#################################################################################
.globModel <- function(params) {
  switch(params$design,
         PEP_INTENSITY = {
           switch(params$residual.variability,
                  biological = {
                    model = "M~1+sample+pep"
                  },
                  technical = {
                    model = "M~1+sample"
                  }
           )
           if(params$channel==1){model=paste0(model,"+channel")}
         } ,
         PEP_INTENSITY = {
           switch(params$residual.variability,
                  biological = {
                    model = "M~1+sample+pep"
                  },
                  technical = {
                    model = "M~1+sample"
                  }
           )
           if(params$channel==1){model=paste0(model,"+channel")}
         } ,
         PEP_RATIO = {
           switch(params$residual.variability,
                  biological = {
                    model = "M~1+pep"
                  },
                  technical = {
                    model = "M~1"
                  }
           )
           if(params$channel==1){model=paste0(model,"+channel")}
         }
  )
  return(model)
}
#################################################################################
## Name : .analysisDiff
## Arguments:
##           "object" data.frame
##           "parameters" data.frame
## Description :
#################################################################################
.analysisDiff <- function(object,params)
{
  print("Execution of .analysisDiff")
  
  if(params$residual.variability == "biological" & length(unique(object$repTech))!=1){ 
    
    print("Technical Replicates are being merged ...")
    objTec = object %>% group_by(proteinId, peptide, experiment, sample,  replicate ) %>% 
      dplyr::summarise(M = mean(M), A = mean(A), 
                       quantifSet = paste(quantifSet,collapse="+"),
                       peptideId = paste(peptideId,collapse="+")) %>% ungroup()
    
    objTec = left_join( objTec, object )
    objTec$repTech <- "RepTech"
    object <- objTec  
  }
  
  if(params$design == "PEP_INTENSITY"){ # Case PEP_INTENSITY
    object = object %>% mutate(sample=as.factor(sample)) %>%  group_by(proteinId)
    .lmRepParam = function( x ){ .lmProteinPEP_INTENSITY( subObject = x , params = params )}
    
  }else if(params$design=="PEP_RATIO"){ # Case PEP_RATIO
    conditions = object$sample %>% unique
    if(length(conditions)==1){
      pairs = as.matrix(conditions)
    }else{
      pairs = conditions %>% combn(2)
    }
    # Build the data for the secondary ratios
    tmpPairs = NULL
    for(i in 1 : dim(pairs)[2])
    {
      tmp0 = paste(pairs[,i],collapse="//")
      tmpPairs = tmpPairs %>% rbind(object %>% filter(sample %in% pairs[,i]) %>% mutate(cond=tmp0) %>% group_by(proteinId) %>% nest %>% mutate(condId=tmp0) ) 
    }
    # Build the data for the primary ratios
    tmpCond = NULL
    for(i in 1 : length(conditions))
    {
      tmpCond = tmpCond %>% rbind( object %>% filter(sample%in%conditions[i]) %>% 
                                     mutate(cond=conditions[i]) %>% group_by(proteinId) %>% 
                                     nest %>% mutate(condId=conditions[i]) 
      )
    }
    object = tmpPairs %>% full_join(tmpCond,b=c("proteinId","condId")) %>% 
      mutate(data.x.NULL=map(data.x,is.null)) %>% 
      unnest(data.x.NULL) %>% mutate(data=data.x)
    object$data[object$data.x.NULL] = object$data.y[object$data.x.NULL]
    object = object %>% select(proteinId,condId,data) %>% unnest() %>% group_by(proteinId,condId)
    
    .lmRepParam = function(x) .lmProteinPEP_RATIO(x,params=params)
  }else{stop("wrong name of design")}
     
  object = object %>% do(.lmRepParam(.)) %>% ungroup()
  if(params$design=="PEP_RATIO"){
    object = object %>% select(-condId)
   }
  
  # # start DEBUG chapter
  # for(i in 1:length(object$data)){ # DEGUG
  #   subObject = object$data[[i]]# DEGUG
  #   subRes = .lmRepParam(subObject) # DEBUG
  #   print( paste( "i = " , i ) ) # DEGUG
  #   print( paste( "Number of peptides : ", length(unique(subObject$peptide) ) ) )
  #   print(".lmRepParam(subObject) = " ) # DEGUG
  #   print( subRes ) # DEGUG
  #   print("####################################") # DEGUG
  #   # # Print withe the correction of the effects
  #   #  eq = as.formula(subRes %>% filter(term=="State2.State1") %>% .$model)
  #   #  effects = eq %>% as.character() %>% .[3] %>% strsplit(.,split=" \\+ ") %>% unlist %>% .[-1]
  #   #  lm0 = lm(eq,subObject)
  #   #  # Build the table of prediction
  #   #  predVSraw = subObject %>% rename(peptideOld=peptide,replicateOld=replicate) %>%
  #   #    mutate(peptide = subObject$peptide %>% unique %>% .[order(.)] %>% .[1],
  #   #           replicate = subObject$replicate %>% unique %>% .[order(.)] %>% .[1] ) %>%
  #   #    mutate(type="predicted") %>% mutate(M=predict(lm0)) %>% mutate(peptide=peptideOld,replicate=replicateOld) %>%
  #   #    select(-replicateOld,-peptideOld) %>% full_join( (subObject %>% mutate(type="raw")))
  #   # 
  #   #  # Plot the raw data and the corrected data thanks to the adjustement found in lm
  #   #  predVSraw %>% unite(truc,c("type","sample"),remove=FALSE) %>%
  #   #    ggplot(aes(peptide,M,color=replicate)) + geom_point() + facet_wrap(~truc)
  # 
  # # if(length(unique(subObject$sample))==1){ # Don't plot the results # DEBUG
  # #   if(parameters$V2[parameters$V1=="design"]=="LABELFREE"){subRatio = subRes[1:dim(pairs)[2],];subRes=subRes[(dim(pairs)[2]+1):dim(subRes)[1],]} # DEGUG
  # #   p=.myBoxPlot(rawDataPep = subObject,ratios =subRes ) # DEGUG
  # #   if(parameters$V2[parameters$V1=="design"]=="LABELFREE"){ # DEGUG
  # #     library(cowplot) # DEGUG
  # #     p1=subRatio %>% ggplot(aes(term,estimate)) + geom_point(aes(color=-log10(p.value),size=-log10(p.value))) + # DEGUG
  # #       geom_errorbar(aes(ymin=ci2.5,ymax=ci97.5),linetype=3,width=0.5) + ylab("log2 ratio") # DEGUG
  # #     print(plot_grid(p,p1)) # DEGUG
  # #   }else{ # DEGUG
  # #     print(p) # DEGUG
  # #   } # DEBUG
  # # } # DEGUG
  # # scan() # DEGUG
  # # if(length(unique(subObject$sample))<2){scan()}
  # }  # DEBUG
  # 
  # # End DEBUG chapter
  
  # Correction of the p.values if needed
  .myAdjust = function(x){p.adjust(x$p.value,method = params$pAdj.method)}
  object = object %>% group_by(term) %>% nest %>%
    mutate(p.valueAdjusted = map(data,.myAdjust )) %>% unnest() %>%
    dplyr::rename(p.valueNonAdjusted=p.value) %>% dplyr::rename(p.value = p.valueAdjusted)   
 
  ind <- which(object$std.error==0) 
  if (length(ind) ==  0 ) {
  
  # Correction of the CI
  .myCorrectionCI <- function(x){
    n = dim(x)[1]
    alpha = 0.05/n
    x = x %>% mutate( ci2.5 = estimate+std.error*qt(1-alpha,df), ci97.5 = estimate+std.error*qt(alpha,df) ) 
    return(x)
  }
  
  if(params$pAdj.method!="none" & length( which( table(as.data.frame(object)$term)==1) ) ==0){
    object = object %>% group_by(term) %>% nest %>% 
      mutate(data = map(data,.myCorrectionCI)) %>% unnest()
  }
  
  # Remove NaN
  numericColumn =c("p.value","estimate","std.error","p.valueNonAdjusted","ci2.5","ci97.5")
  for(i in numericColumn){
    object[is.nan(unlist(object[,i])),i] = NA
    }
  }  
  resu = object 
  print(".analysisDiff done")
  return(resu) 
}


#################################################################################
## Name : .numOrder
## Arguments:
## Description :
#################################################################################
.numOrder <- function(x){
  if(length(x)==1){return(1)}
  gr <- gregexpr( "[0-9\\.]+" , x ) 
  gr <- sapply( regmatches( x , gr ) , as.numeric ) %>% t
  perm = 1:dim(gr)[1]
  for( i in dim(gr)[2]:1 ){
    permTmp = order( gr[,i] , decreasing = FALSE )
    gr = gr[permTmp,]
    perm = perm[ permTmp ]
  }
  return(perm)
}

#################################################################################
## Name : .grepNum
## Arguments:
## Description :
#################################################################################
.grepNum <- function(x){
  gr <- gregexpr("[0-9\\.]+" , x )
  gr <- sapply(regmatches(x , gr) , as.numeric) %>% t %>% 
    apply(MARGIN = 1 ,FUN= function(x) paste( x , collapse = "" ) )
  return(gr)
}

#################################################################################
## Name : .lmProteinPEP_INTENSITY
## Arguments:
##           "subObject" data.frame
##           "residualVariability" data.frame
## Description :
#################################################################################
.lmProteinPEP_INTENSITY<-function(subObject,params){   
  subObject = subObject %>% arrange(sample)  
  # Compute the effectiv of each possible terms of the equantion for the linear model
  vect = apply(subObject[,!names(subObject) %in% c(params$residual.variability,"A","M")],MARGIN=2,function(x)length(unique(x))) %>% 
    t %>% as_tibble
  conditions = subObject$sample %>% unique
  conditions = conditions[ sapply(conditions,.grepNum) %>% as.numeric() %>% order(decreasing = TRUE)] %>% as.character
  
  # Select the correct fix effect
  if(params$residual.variability=="biological"){
    fixEffect = c("sample","experiment","peptide") # ,"replicate"
  }else if(params$residual.variability=="technical"){
    fixEffect = c("sample","experiment")
  }else{ stop("wrong name of design") }
  
  # Equation
  eq = "M~1"
  for(i in fixEffect){
    eq = paste0(eq,ifelse(vect[,i]>1,paste0("+",i,sep=""),""))
  }
  
  # Model primary ratio
  if( subObject$sample %>% unique %>% length() == 1 ){
    model = lm(as.formula(eq),subObject)
    CI = model %>% confint ; CI = CI %>% as_tibble %>% mutate( term = rownames(CI) ) ; colnames(CI) = c("ci2.5","ci97.5","term")
    outputPrim = model  %>% tidy %>% full_join(CI,by="term") %>% filter(term=="(Intercept)") %>% 
      mutate( term = unique(subObject$sample) , model = eq , term = as.character(term) )
    outputPrim$df = summary(model)$df[2]
  }else{
    # Apply a linear model with the model defined by "eq" on the protein, then extract the mean value of State1 & State2 
    eqPrimary = gsub(pattern = "1",replacement = "-1",x = eq)
    model = lm(as.formula(eqPrimary),subObject)
    CI = model %>% confint ; CI = CI %>% as_tibble %>% mutate( term = rownames(CI) ) ; colnames(CI) = c("ci2.5","ci97.5","term")  
    outputPrim = model  %>% tidy %>% full_join(CI,by="term") %>% mutate( keep =  str_detect(term,pattern = "sample") ) %>% 
      filter(keep) %>% select(-keep) %>% mutate(term = map(term,function(x) gsub("sample","",x) )) %>% unnest(term) %>% 
      mutate( term = as.character(unique(subObject$sample)) , model = eqPrimary  )
    outputPrim$df = summary(model)$df[2]
    # Put NA on p.value, std.error and ci if there is only one replicate observed in a sample
    if(params$residual.variability=="biological"){
      outputPrim = subObject %>% group_by(sample) %>% nest %>% mutate(nRep=map(data,function(x)length(unique(x$replicate))) ) %>%
        select(-data) %>% unnest() %>% dplyr::rename(term=sample) %>% full_join(.,outputPrim,by="term") %>% 
        mutate(p.value=ifelse(nRep<=1,NA,p.value),
               std.error=ifelse(nRep<=1,NA,std.error),
               ci2.5=ifelse(nRep<=1,NA,ci2.5),
               ci97.5=ifelse(nRep<=1,NA,ci97.5)
               ) %>% select(-nRep)
    }
  }
  
  # Model secondary ratio
  outputSec = tibble(ref = conditions[1:(length(conditions)-1)],data=conditions[1:(length(conditions)-1)]) %>% group_by(ref) %>% nest
  if(subObject$sample %>% unique %>% length() ==1){
    outputSec = tibble(term="remove",estimate=1,std.error=1,statistic=1,p.value=1,model="1",ci2.5=1,ci97.5=1,df=1)
  }else{
    .tmp1LmProteinPEP_INTENSITY <- function(x) .lmByRef(x,eq,subObject,params)
    outputSec = outputSec %>% mutate( data = map( data , .tmp1LmProteinPEP_INTENSITY ) ) %>% unnest(data) %>% select(-ref)
  }
  
  # Merge the 2 outputs
  output = outputPrim %>% full_join(outputSec,by=c("term","estimate","std.error","statistic","p.value","ci2.5","ci97.5","model","df"))

  # Fill the NA and the Inf
  naAndInf = subObject$sample %>% levels %>% .[order(as.numeric(sapply(.,.grepNum)),decreasing = TRUE)] %>% 
    combn(2) %>% t %>% as_tibble %>% unite(term,V1,V2,sep=".") %>% 
    filter(!term%in%output$term) 
  # When there is no NA of Inf
  if( dim( naAndInf )[1]!=0 ){ 
    .tmp2LmProteinPEP_INTENSITY <- function(x) .fillNaPEP_INTENSITY(x,subObject)
    naAndInf = naAndInf %>% mutate( data = map( term , .tmp2LmProteinPEP_INTENSITY ) ) %>% unnest(data) %>% select(-term1) %>%
      mutate_if(is.logical,as.double)
    naAndInf$df = NA
  }else{
    naAndInf = tibble(term="remove",estimate=1,std.error=1,statistic=1,p.value=1,model="1",ci2.5=1,ci97.5=1,df=0) 
  }
  
  # Merge the results
  output = output %>% full_join(naAndInf,by=c("term","estimate","std.error","statistic","p.value","model","ci2.5","ci97.5","df")) %>% filter(term!="remove") %>% 
    mutate(df=as.numeric(df))
  return(output)
}
#################################################################################
## Name : .fillNaPEP_INTENSITY
## Arguments:
## Description :
################################################################################
.fillNaPEP_INTENSITY <- function(term,subObject){
  term1 = strsplit(term,split="\\.") %>% unlist %>% .[1]
  term2 = strsplit(term,split="\\.") %>% unlist %>% .[2]
  m1 = subObject %>% filter(sample==term1) %>% .$M %>% mean() %>% ifelse(is.nan(.),0,.)
  m2 = subObject %>% filter(sample==term2) %>% .$M %>% mean() %>% ifelse(is.nan(.),0,.)
  estimate = log2(m1/m2) %>%  ifelse(is.nan(.), NA,.)
  output = tibble(term=term,estimate=estimate,std.error=NA,statistic=NA,p.value=NA,ci2.5=NA,ci97.5=NA) %>% 
    mutate( std.error = as.numeric( std.error ),
            statistic = as.numeric( statistic ),
            p.value = as.numeric( p.value ),
            ci2.5 = as.numeric(ci2.5),
            ci97.5 = as.numeric(ci97.5),
            model="specific" )
  return(output)
}

#################################################################################
## Name : .lmByRef
## Arguments:
## Description :
################################################################################
.lmByRef = function(ref,eq,subObject,params){
  # Compute the model
  ref = as.character(ref)
  iN = as.numeric(.grepNum(ref))
  subObject$sample = relevel(as.factor(subObject$sample), ref=ref)
  model = lm(as.formula(eq),subObject)
  tmp = model %>% tidy
  
  # Get the parameters of the model
  tmp = tmp[grep("sample",tmp$term),] 
  tmp = tmp %>% mutate( numSample = map( term , .grepNum ) ) %>% unnest(numSample) %>% 
    mutate(numSample=as.numeric(numSample))
  
  # Confint
  cI = model %>% confint
  rNCI = rownames(cI)
  cI = cI %>% as_tibble %>% mutate(term=rNCI) 
  colnames(cI)=c("ci2.5","ci97.5","term")
  
  tmp = tmp %>% full_join(cI,by="term")
  
  # Build the output
  output = tmp %>% filter(numSample<iN) %>% select(term,estimate,std.error,statistic,p.value,ci2.5,ci97.5)
  output$term = gsub(pattern="sample",paste0(ref,"."),output$term)
  output$model = eq
  output$df = summary(model)$df[2]
  
  # Do not report p.value or variance if the number of replicates obeserved in once condition is <1
  if(params$residual.variability=="biological"){
    nRep = subObject %>% group_by(sample) %>% nest %>% mutate(nRep=map(data,function(x)length(unique(x$replicate))) ) %>%
      select(-data) %>% unnest() %>% .$nRep
    if(!prod(nRep>1)){
      output$p.value=NA
      output$std.error=NA
      output$ci2.5=NA
      output$ci97.5=NA
    }
  }

  
  # Invert the sign of S2.S1 because the model compute S1.S2
  output[,c("estimate","ci2.5","ci97.5")] = -1* output[,c("estimate","ci97.5","ci2.5")]
  return(output)
}
#################################################################################
## Name : .primaryRatiosPEP_RATIO
## Arguments:
## Description :
#################################################################################
.primaryRatiosPEP_RATIO<- function(subObject,params){
  # Compute the effectiv of each possible terms of the equantion for the linear model
  vect = apply(subObject[,!names(subObject) %in% c(params$residual.variability,"A","M")],MARGIN=2,function(x)length(unique(x))) %>% 
    t %>% as_tibble
  
  # Select the correct fix effect
  if(params$residual.variability=="biological"){
    fixEffect = c("experiment","peptide")
  }else if(params$residual.variability=="technical"){
    fixEffect = NULL
  }else{ stop("wrong name of design") }
  
  # Equation
  eq = "M~1"
  for(i in fixEffect){
    eq = paste0(eq,ifelse(vect[,i]>1,paste0("+",i,sep=""),""))
  }
  
  # Model
  myLM = lm(as.formula(eq),subObject)
  
  # Output
  ci = myLM %>% confint() %>% .[1,] %>% t %>%  as_tibble()
  names(ci) = c("ci2.5","ci97.5")
  coef = myLM %>% tidy 
  output = cbind(coef[1,],ci) %>% as_tibble
  output$term = unique(subObject$sample)
  output$model = eq
  output$df = myLM %>% summary %>% .$df %>% .[2]
  return(output)
} 
#################################################################################
## Name : .secondaryRatiosPEP_RATIO
## Arguments:
## Description :
#################################################################################
.secondaryRatiosPEP_RATIO <- function(subObject,params){
  # Compute the effectiv of each possible terms of the equantion for the linear model
  vect = apply(subObject[,!names(subObject) %in% c(params$residual.variability,"A","M")],MARGIN=2,function(x)length(unique(x))) %>% 
    t %>% as_tibble
  conditions = subObject$cond %>% unique %>% strsplit(split="//") %>% unlist
  
  if(unique(vect$sample)==1){
    conditions = conditions[conditions %>% .grepNum() %>% order(decreasing = TRUE)]
    m1 = subObject %>% filter(sample==conditions[1]) %>% .$M %>% mean %>% ifelse(is.na(.),0,.)
    m2 = subObject %>% filter(sample==conditions[2]) %>% .$M %>% mean %>% ifelse(is.na(.),0,.)
    output = tibble(
      term = paste( conditions , collapse = "." ),
      estimate = log( m1 / m2 ),
      std.error = NA,
      statistic = NA,
      p.value = NA,
      ci2.5 = NA,
      ci97.5 = NA,
      model = "specific",
      df = 0
    )
    return(output)
  }
  
  # Select the correct fix effect
  if(params$residual.variability=="biological"){
    fixEffect = c("experiment","peptide")
  }else if(params$residual.variability=="technical"){
    fixEffect = "experiment"
  }else{ stop("wrong name of design") }
  
  # Equation
  eq = "M~1+sample"
  for(i in fixEffect){
    eq = paste0(eq,ifelse(vect[,i]>1,paste0("+",i,sep=""),""))
  }
  
  # Model
  myLM = lm(as.formula(eq),subObject)
  
  
  # Output
  coef = myLM %>% tidy 
  cI = myLM %>% confint %>% .[1,] 
  cI = cI %>% t %>%  as_tibble()
  names(cI) = c("ci2.5","ci97.5")
  
  
  output = cbind(coef[1,],cI) %>% as_tibble
  output$term = paste(conditions[conditions %>% .grepNum %>% order(decreasing = TRUE)],collapse=".")
  output$model = eq
  output$df = myLM %>% summary %>% .$df %>% .[2]
  
  # Put NA in p.value and variance if a condition have been observed only in one run
  if(params$residual.variability=="biological"){
    nRep = subObject %>% group_by(sample) %>% nest %>% mutate(nRep=map(data,function(x) x %>% .$replicate %>% unique %>% length)) %>% 
      select(-data) %>% unnest() %>% .$nRep
    if(!prod(nRep>2)){
      output$std.error=NA
      output$statistic=NA
      output$p.value=NA
      output$ci2.5=NA
      output$ci97.5=NA
      output$df=NA
    }
  }
  
  # Invert the estimate because we compute the estimate in the wrong way 
  output[,c("estimate","ci2.5","ci97.5")] = -1* output[,c("estimate","ci97.5","ci2.5")]
  return(output)
} 
#################################################################################
## Name : .lmProteinPEP_RATIO
## Arguments:
##           "subObject" data.frame
##           "residualVariability" data.frame
## Description :
################################################################################
.lmProteinPEP_RATIO<-function(subObject,params){
  n = subObject$cond %>% unique %>% strsplit(split="//") %>% unlist %>% length
  if( n == 1 ){ output = .primaryRatiosPEP_RATIO( subObject , params ) }
  if( n == 2 ){ output = .secondaryRatiosPEP_RATIO( subObject , params )}
  return(output)
}

#################################################################################
## Name : .MAToIntensity
## Arguments:
##           "object" data
## Description :
#################################################################################
.MAToIntensity <- function(object){
  # Rebuild the initial data shape of the nonNorm data
  tmp = object %>% mutate(intensityL = (M+2*A)/2,intensityR = (2*A-M)/2) %>%
    select(proteinId,peptide,experiment,replicate,repTech,sample,-M,-A,intensityL,intensityR) %>%
    separate(sample,c("sampleL","sampleR"),sep="\\.")
  tmp = rbind( select(tmp,-sampleL,-intensityL) %>% dplyr::rename(intensity=intensityR) %>% dplyr::rename(sample=sampleR) ,
               select(tmp,-sampleR,-intensityR) %>% dplyr::rename(intensity=intensityL) %>% dplyr::rename(sample=sampleL) ) %>%
    group_by(proteinId,peptide,experiment,replicate,repTech,sample) %>% dplyr::summarise(intensity=mean(intensity)) %>%
    ungroup
  return(tmp)
}

#################################################################################
## Name : .cvOutliers
## Arguments:
##           "X" list of 2 tables :
##                   "data" the data without the outliers dected previously
##                   "outlier" the list of outliers detected previously, if no
##                   outliers have been detected add NULL.
##           "rep" if "bio" remove the biological outliers, if "tech" remove the
##                 technical outliers
## Description : Remove the biological or thechnical replicate based on the CV.
#################################################################################
.cvOutliers = function(X,rep)
{
  switch(rep,
         bio = {
           rep="repTech" # Not a mistake : group_by(repTech) then summarise to compute the CV
           outName="outBioRep" # The name of the outlier for tracability
         },
         tech = {
           rep="replicate" # Not a mistake : group_by(replicate) then summarise to compute the CV
           outName="outTechRep" # The name of the outlier for tracability
         }
  )
  # Find the outlier
  tmp = X$data %>%
    group_by_(.dots=c("sample","experiment","proteinId","peptide",rep)) %>% dplyr::summarise(CV=sd(M)/mean(M)) %>% ungroup
  # Compute the CV on the bio or tech replicate
  tmp$out = tmp$CV %in% boxplot(tmp$CV,plot=FALSE)$out # Find the outliers
  # Table of the outliers
  outlier = tmp %>% filter(out==TRUE) 
  if(dim(outlier)[1]!=0){
    outlier = outlier %>% mutate( out = "outBioRep" ) %>% # Get the outliers
      left_join(X$data,by=c("sample","experiment","proteinId","peptide",rep)) %>% select(-CV) # join to recover the M & A columns
  }
  
  # Table of the data
  data =  tmp %>% filter(out==FALSE) %>% select(-CV,-out) %>% # Remove the outliers from the data
    left_join(X$data,by=c("sample","experiment","proteinId","peptide",rep)) # Join to recover the M & A columns
  # Fill the outlier table in X if NULL
  if(is.null(X$outlier)){
    X$outlier=outlier[0,]
  }
  # Join the old table of the outliers with the new one.
  outlier = rbind(outlier,X$outlier)
  # Output
  X = list(
    data = data ,
    outlier = outlier )
  return(X)
}

#################################################################################
## Name : .peptideOutliers
## Arguments:
##           "X" list of 2 tables :
##                   "data" the data without the outliers dected previously
##                   "outlier" the list of outliers detected previously, if no
##                   outliers have been detected add NULL.
## Description : Remove the peptides outliers in each proteins.
#################################################################################
.peptideOutliers = function(X)
{
  .fct <- function(x){ x$M %in% boxplot(x$M,plot=FALSE)$out }
  # Find the outlier
  tmp = X$data %>%
    group_by(sample,experiment,proteinId,repTech,replicate) %>% nest %>% mutate(out = map(data,.fct) ) %>% unnest()
  # Table of the outliers
  outlier = tmp %>% filter(out==TRUE)
  if(dim(outlier)[1]!=0){
    outlier = outlier %>% mutate( out = "outPeptide" )
  }
  # Table of the data
  data =  tmp %>% filter(out==FALSE) %>% select(-out)
  # Fill the outlier table in X if NULL
  if(is.null(X$outlier)){
    X$outlier=outlier[0,]
  }
  # Join the old table of the outliers with the new one.
  outlier = rbind(outlier,X$outlier)
  # Output
  X = list(
    data = data ,
    outlier = outlier )
  return(X)
}

#################################################################################
## Name : .mixBioPeptideOutliers
## Arguments:
##           "X" list of 2 tables :
##                   "data" the data without the outliers dected previously
##                   "outlier" the list of outliers detected previously, if no
##                   outliers have been detected add NULL.
## Description : Remove the peptides outliers in each proteins.
#################################################################################
.mixBioPeptideOutliers <- function(X,params){
  .fct <- function(x){ x$M %in% boxplot(x$M,plot=FALSE)$out }
  # Normalize the value of each observation in the protein by the median of the 
  # observations across the replicates (tech or bio depending on the parameters)
  rep = ifelse(params$residual.variability=="biological","repTech","replicate") # not a bug, the median is not computed on rep
  tmp = X$data %>%
    group_by_( .dots = c("sample","experiment","proteinId","peptide",rep) ) %>% dplyr::summarise(m = median(M,na.rm=TRUE))
  tmp = X$data %>% full_join(tmp,by=c("sample","experiment","proteinId","peptide",rep)) %>% mutate(M=M/m) %>% select(-m) %>% 
    group_by(sample,experiment,proteinId) %>% nest()  %>% mutate(out = map(data,.fct) ) %>% unnest() #unnest(cols = c(data, out))
  
  # Table of the outliers : remove the extreme values thanks to a boxplot method.
  outlier = tmp %>% filter(out==TRUE)
  if(dim(outlier)[1]!=0){ # If outliers are detected by the method of the boxplot, get the values from the data, join two table, label the outliers and continue searching for others outliers
    outlier = outlier %>% select(-M,-out) %>% left_join((X$data),by=names(.)) %>% mutate( out = "outMixBioPeptide" ) 
  }
  # Table of the data
  data =  tmp %>% filter(out==FALSE) %>% select(-out,-M) %>% left_join(X$data,by=names(.))
  
  # Force every peptides to be present in both states, else flag the peptide as an outlier
  dataPep = data %>% select(sample,experiment,proteinId,peptide) %>% unique %>% mutate(out="inData") %>% spread(sample,out) %>% 
    group_by(experiment,proteinId,peptide) %>% nest()   %>% 
    mutate(out=map(data,function(x) ifelse( (x %>% as.matrix %>% .[!is.na(.)] %>% length)<2 ,"outBalanceStates","FALSE") )) %>% 
    select(-data) %>% unnest() #unnest(cols = c(out))
  dataOutput = data %>% full_join(dataPep) %>% filter(out =="FALSE") %>% select(-out)
  # if there is no outliers do not join the new outlier table with the old one, if not do the join
  if(dim(outlier)[1]==0){
    outlier = data %>% full_join(dataPep) %>% filter(out == "outBalanceStates") 
  }else{
    outlier = data %>% full_join(dataPep) %>% filter(out == "outBalanceStates") %>% full_join(outlier)
  }
  data = dataOutput
  
  # Fill the outlier table in X if NULL
  if(is.null(X$outlier)){
    X$outlier=outlier[0,]
  }
  # Join the old table of the outliers with the new one.
  outlier = bind_rows(outlier,X$outlier)  # outlier = rbind(outlier,X$outlier)
  # Output
  X = list(
    data = data ,
    outlier = outlier )
  return(X)
}

#################################################################################
## Name : .mySpread
## Arguments:
## Description :
#################################################################################
.mySpread = function(data,value,newName,toSpread="term",key="ProteinID",sep=""){
  tmp = data %>% select_(key,toSpread,value) %>% dplyr::rename_(.dots = setNames(toSpread, newName)) %>%   spread_(newName,value,sep=sep)
  return(tmp)
}

#################################################################################
## Name : .myBoxPlot
## Arguments:
## Description :
#################################################################################
.myBoxPlot <- function( rawDataPep ,ratios ){
  p = rawDataPep %>% ggplot( aes( sample, M ) ) +
    geom_boxplot( aes( fill = sample ) ) +
    geom_point( position = position_dodge ( width = 0.75 ), aes( group = sample , color = peptide) ) + scale_fill_brewer(palette="Dark2") +
    geom_point( data = ratios, aes( term, estimate,size = -log10(p.value) ) ) + labs(y="log2 intensity",x="condition")
  return(p)
}

#################################################################################
## Name : .control
## Arguments: object
##            params
##            normProt
## Description : check the frequent errors in the data, the parameters or in the
##               protein of normalization.
#################################################################################
.control <- function(object,params,normProt,dataRef){
  
  message = NULL
  # ---------- Control the correct number and type of files
  condition = sum(exists("object"), exists("params"))
  if( !condition ){
    message = "files table.txt or param_char.txt (or both) is(are) missing"
    stop(message)
  }
  
  condition = sum(exists("dataRef"),c("normalization.ref.test") %in% names(params))
  if( prod(condition, c("normalization.ref.test") %in% names(params)) == 1){ 
    message = "file tableRef.txt is missing"
    stop(message)
  }
  
  # ---------- Control the correct type of the parameters
  condition = prod(c("design","normalization.method") %in% names(params))
  if( !condition ){
    message = "wrong number or names of parameters, parameters must contains:\n              design\n              normalization.method\n              pAdj.method\n              residual.variability"
    stop(message)
  }
  
  # design
  p = c("PEP_INTENSITY","PEP_RATIO")
  condition = params$design %in% p
  if( !condition ){
    message = paste(c("wrong name of design,\n         parameters allowed :",p),collapse="\n                             ")
    stop(message)
  }
  
  # normalization.method
  p =  c("none.none",
         "none.scale",
         "median.none",
         "median.scale",
         "mean.none",
         "mean.scale",
         "quantile")
  condition = params$normalization.method %in% p
  if( !condition ){
    message = paste(c("wrong name of normalization.method,\n         parameters allowed :",p),collapse="\n                             ")
    stop(message)
  }
  
  # normalization.ref.test
  p =  c("none.none.mean.var",
         "none.none.mean.scale",
         "none.none.median.var",
         "none.none.median.scale",
         "none.scale.mean.var" ,
         "none.scale.mean.scale",
         "none.scale.median.var",
         "none.scale.median.scale",
         "median.none.mean.var",
         "median.none.mean.scale",
         "median.none.median.var",
         "median.none.median.scale",
         "median.scale.mean.var",
         "median.scale.mean.scale",
         "median.scale.median.var",
         "median.scale.median.scale",
         "mean.none.mean.var",
         "mean.none.mean.scale",
         "mean.none.median.var",
         "mean.none.median.scale",
         "mean.scale.mean.var",
         "mean.scale.mean.scale",
         "mean.scale.median.var",
         "mean.scale.median.scale",
         "quantile.none.mean.var",
         "quantile.none.mean.scale",
         "quantile.none.median.var",
         "quantile.none.median.scale"  )
  if ("normalization.ref.test" %in% names(params)) condition = params$normalization.ref.test %in% p
  if ( !condition ){
    message = paste(c("wrong name of normalization.ref.test,\n         parameters allowed :",p),collapse="\n                             ")
    stop(message)
  }
  
  # pAdj.method
  if ("pAdj.method" %in% names(params)) { 
            p =  p.adjust.methods
            condition = params$pAdj.method %in% p}
  if( !condition ){
    message = paste(c("wrong name of pAdj.method,\n         parameters allowed :",p),collapse="\n                             ")
    stop(message)
  }
  
  # residual.variability
  if ("residual.variability" %in% names(params)) {
             p =  c("biological","technical")
             condition = params$residual.variability %in% p }
  if( !condition ){
    message = paste(c("wrong name of residual.variability,\n         parameters allowed :",p),collapse="\n                             ")
    stop(message)
  }
  
  # ---------- Control the correct type of the data
  # Check the fullness of the data
  condition = object %>% filter(!is.na(value)) %>% dim %>% .[1]
  if( !condition ){
    message = "no intensity values are observed in table.txt."
    stop(message)
  }
  
  # Check the number of conditions
  condition = (object %>% filter(!is.na(value)) %>% .$sample %>% unique %>% length )-1
  if( !condition ){
    message = "onely one condition is observed in table.txt."
    stop(message)
  }
  
  # ---------- Control the correct type of the normalization protein
  condition = sum(normProt$proteinId %in% (object %>% .$proteinId %>% unique)) || is.null(normProt)
  if( !condition ){
    message = "wrong names of the protein of normalization."
    stop(message)
  }
  
  
  print("Control done")
  return(1)
}

#################################################################################
## Name : normEachProt
## Arguments: object
##            objectRef
##            params
## Description : normalize each protein by the mean & the variance of the twin 
##               protein (with the same proteinId)  found in objectRef. The 
##               parameters of normalization are fount in params.
#################################################################################

.normEachProt <- function(object,objectRef,params){
   
  # dispersion estimated by quantile 
  quant <-function(REF,DAT){
    bx <-boxplot(REF,plot=FALSE)
    if (length(bx$out)!=0){
      max <- bx$stats[5]
      min <- bx$stats[1]
      REF <- REF[ REF > min & REF < max ]
    }
    sdx <-sd(REF, na.rm = FALSE)
    sdy <-sd(DAT, na.rm = FALSE)
    X.to.determine.normalization <- matrix(REF)
    Y.to.normalize               <- matrix(DAT)
    
    if ( is.na(sdx)  | is.na(sdy) ){ 
      Y.normalized <- Y.to.normalize 
    } else if ( sdx <= sdy ) {       
      Y.normalized <- Y.to.normalize               
    } else if ( sdx > sdy ) {        
      target <- normalize.quantiles.determine.target(X.to.determine.normalization,
                                                     dim(X.to.determine.normalization)[1])   
      Y.normalized <- normalize.quantiles.use.target(Y.to.normalize,target) 
    } 
    return(Y.normalized)
  }  #quant <-function(REF,DAT){

# effect 
.effectByProtINT <- function(object,objectRef){

objectNorm <-NULL   

for (i in 1:length(unique(object$proteinId))){ 
  prot         <- unique(object$proteinId)[i]
  ind          <- which(object$proteinId == prot)
  subObject    = object[ind,]
  subObjectRef = objectRef %>% filter(proteinId == 
                                        unlist(str_split(unique(subObject$proteinId), "\\-"))[1])
  
  namedat         <- names(table(subObjectRef$sample))
  nameref         <- names(table(   subObject$sample))
  namediff1     <-setdiff(unique(c(  namedat, nameref)),intersect( namedat, nameref) ) # les states qui ne matchtent pas 
  namediff2     <-setdiff(unique(c(  namedat, nameref)),namediff1)

  matrixdat       <- as.data.frame(subObject)      
  matrixref       <- as.data.frame(subObjectRef)   
  matrix          <- data.frame( type   = c(rep("data",   nrow(matrixdat)), 
                                            rep("dataref",nrow(matrixref))    ),
                                 sample = c(as.character(      matrixdat$sample),
                                            as.character(      matrixref$sample) ),
                                 peptide= c(as.character(      matrixdat$peptide),
                                            as.character(      matrixref$peptide)),
                                 M      = c(                   matrixdat$M,
                                                               matrixref$M)     )
  
  sds           <- aggregate(M ~ paste(type,sample,sep="_") , matrix, FUN = function(x) {sd(x,na.rm=TRUE)})
  colnames(sds) <- c("type*state","sd")
  sds$sample    <- gsub(".*_","",sds$'type*state')
  S <-NULL
  for (s in namediff2){
    S <-c(S,prod(sds$sd[which(sds$sample==s)]))
  }
  inddif        <- c(namediff1, namediff2[which(S==0)] )
  namediff3     <-setdiff(namediff2,inddif)
  lendif <- length(namediff3)
  lenuni <- length(unique(c(  namedat, nameref)))
  
  if (dim(subObjectRef)[1] ==0  |                                               # si la proteine n'existe pas dans dataref
      sum(is.element("State1",nameref), is.element("State1",namedat)   )!=2  |  # si le State1 n'est pas dans data et dataref
      prod(sds$M[which(sds$sample=="State1")]  )==0                          |  # si le State1 ne varie pas dans data et dataref
      lendif<2                                                               |
      prod(is.nan(subObjectRef$M))==1                                        |
      prod(is.na( subObjectRef$M))==1                                        )
  {
    subObjectNorm   <- subObject
    subObjectNorm$M <- rep(NA,dim(subObject)[1])
    
  } else {   
    
    if (lendif > 2 & lendif < lenuni){ 
      subObject     <- subObject    [-which(is.element(   subObject$sample,inddif)),]
      subObjectRef  <- subObjectRef [-which(is.element(subObjectRef$sample,inddif)),]
     } # if (lendif > 2 & lendif < lenuni)
    
    matrixdat       <- as.data.frame(subObject)      
    matrixref       <- as.data.frame(subObjectRef)   
    matrix          <- data.frame( type   = c(rep("data",   nrow(matrixdat)), 
                                              rep("dataref",nrow(matrixref))    ),
                                   sample = c(as.character(      matrixdat$sample),
                                              as.character(      matrixref$sample) ),
                                   peptide= c(as.character(      matrixdat$peptide),
                                              as.character(      matrixref$peptide)),
                                   M      = c(                   matrixdat$M,
                                                                 matrixref$M)     )
    mean_peptide  = matrix %>%  group_by(type,sample,peptide) %>%  dplyr::summarise(mean.pep = mean(M) ) 
    mean_state    = matrix %>%  group_by(type,sample)         %>%  dplyr::summarise(mean.all = mean(M) )
    mean_peptide  = mean_peptide %>% mutate (type_sample_pep = paste(type,sample,peptide)) %>% 
      mutate (type_sample     = paste(type,sample))
    mean_state    = mean_state   %>% mutate (type_sample     = paste(type,sample))
    matrix_mod   = matrix %>%  mutate( type_sample_pep = paste(type,sample,peptide)) %>% 
      inner_join( mean_peptide, by="type_sample_pep" )     %>% 
      inner_join( mean_state,   by="type_sample" )
    matrix_mod  = matrix_mod %>%  mutate(mean.mod = as.numeric(matrix_mod$mean.pep)  - 
                                           as.numeric(matrix_mod$mean.all) )
    matrix_mod  = matrix_mod %>%  mutate(Mmod = M - mean.mod) 
    matrixmod   = matrix_mod %>%  filter(type.x=="dataref")
    
    ########################          
    # quantile and translation vers les moyennes originales     
    matrixint         <-       matrixdat
    states <- unique(matrix$sample)
    for (s in states){
      inds              <- which(matrixdat$sample==s)
      matrixint$M[inds] <- quant( matrixmod$Mmod[which(matrixmod$sample.x==s)],
                                  matrixdat$M   [which(matrixdat$sample  ==s)])
    }
    Mean1s<-NULL
    for (s in states){
      mean1s            <- mean( matrixint$M[which(matrixint$sample==s)]) - 
        mean( matrixdat$M[which(matrixdat$sample==s)])
      Mean1s <- c(Mean1s,mean1s)
    }   
    MEAN1s <- data.frame(sample=states, mean.sample=Mean1s)
    
    ######################## new        
    matrixapr         <- matrixint
    matrixapr$index <- 1:nrow(matrixapr)          
    mn2       <- merge(matrixapr,MEAN1s, by="sample")
    mn2$M     <- mn2$M - mn2$mean.sample
    matrixapr <- mn2
    matrixapr <- matrixapr[matrixapr$index,]
    nam       <- colnames(matrixint)
    matrixapr <-matrixapr[,nam]

    ######################## new  
    #linear model type for protein effect and sample for condition effect
    model1        = lm(as.formula("M ~ sample*type"),data=matrix)
    matrixeff <- matrixapr
    
    states <- unique(matrixdat$sample)
    for (s in 2:length(states)){
      inds               <- which(matrixapr$sample==states[s])
      matrixeff$M [inds] <-  
      matrixapr$M [inds] - (summary(model1)$coefficients[s]+summary(model1)$coefficients[length(states)+s])   
    }
    subObjectNorm <- matrixeff
    
  } #if (dim(subObject
  objectNorm    <- rbind(objectNorm,subObjectNorm)

} #for (i in )
  
  return(objectNorm)
}  

.effectByProtRAT <- function(object,objectRef){

objectNorm <-NULL   

 for (i in 1:length(unique(object$proteinId))){ 
      prot        <- unique(object$proteinId)[i]
      ind         <- which(object$proteinId == prot)
      subObject   = object[ind,]
      subObjectRef = objectRef %>% filter(proteinId == 
                                            unlist(str_split(unique(subObject$proteinId), "\\-"))[1])
 
      namedat         <- names(table(subObjectRef$sample))
      nameref         <- names(table(   subObject$sample)) 
      namediff1     <-setdiff(unique(c(  namedat, nameref)),intersect( namedat, nameref) ) # les states qui ne machtent pas 
      namediff2     <-setdiff(unique(c(  namedat, nameref)),namediff1)
      matrixdat       <- as.data.frame(subObject)      
      matrixref       <- as.data.frame(subObjectRef)   
      matrix          <- data.frame( type   = c(rep("data",   nrow(matrixdat)), 
                                                rep("dataref",nrow(matrixref))    ),
                                    sample = c(as.character(      matrixdat$sample),
                                               as.character(      matrixref$sample) ),
                                    M      = c(                   matrixdat$M,
                                                                  matrixref$M)     )
      
      sds           <- aggregate(M ~ paste(type,sample,sep="_") , matrix, FUN = function(x) {sd(x,na.rm=TRUE)})
      colnames(sds) <- c("type*state","sd")
      sds$sample    <- gsub(".*_","",sds$'type*state')
  
     S <-NULL
     for (s in namediff2){
      S <-c(S,prod(sds$sd[which(sds$sample==s)]))
     }
    inddif        <- c(namediff1, namediff2[which(S==0)] )
    namediff3     <-setdiff(namediff2,inddif)
    lendif <- length(namediff3)
    lenuni <- length(unique(c(  namedat, nameref)))

  
if (dim(subObjectRef)[1] ==0  |                                                  # si la proteine n'existe pas dans dataref
    #  sum(is.element("State1",nameref), is.element("State1",namedat)   )!=2  |  # si le State1 n'est pas dans data et dataref
    #  prod(sds$M[which(sds$sample=="State1")]  )==0                          |  # si le State1 ne varie pas dans data et dataref
      #lendif<2                                                               |
      prod(is.nan(subObjectRef$M))==1                                        |
      prod(is.na( subObjectRef$M))==1                                        
    )
   {
    subObjectNorm   <- subObject
    subObjectNorm$M <- rep(NA,dim(subObject)[1])
   }  
  else {   
    
    matrixdat       <- as.data.frame(subObject)      
    matrixref       <- as.data.frame(subObjectRef)   
    matrix          <- data.frame( type   = c(rep("data",   nrow(matrixdat)), 
                                              rep("dataref",nrow(matrixref))    ),
                                   sample = c(as.character(      matrixdat$sample),
                                              as.character(      matrixref$sample) ),
                                   M      = c(                   matrixdat$M,
                                                                 matrixref$M)     )

    # quantile and translation vers les moyennes originales     
    ########################
    matrixint         <-       matrixdat
    
    states <- unique(matrix$sample)
    for (s in states){
      inds              <- which(matrixdat$sample==s)
      matrixint$M[inds] <- quant( matrixref$M[which(matrixref$sample.x==s)],
                                  matrixdat$M[which(matrixdat$sample  ==s)])
    }
    
    Mean1s<-Mean2s<-NULL
    for (s in states){
      mean1s            <- mean( matrixint$M[which(matrixint$sample==s)]) - mean( matrixdat$M[which(matrixdat$sample==s)])
      mean2s            <- mean( matrixref$M[which(matrixref$sample==s)]) 
   
      Mean1s <- c(Mean1s,mean1s)
      Mean2s <- c(Mean2s,mean2s)
    }
    MEAN1s <- data.frame(sample=states, mean.sample=Mean1s)
    MEAN2s <- data.frame(sample=states, mean.sample=Mean2s)
    
      #   matrixapr         <- matrixint
      #   matrixapr$M       <- matrixint$M - mean1i 
      #   matrixeff         <- matrixapr
      #   matrixeff$M       <- matrixapr$M - mean2i
      #   subObjectNorm <- matrixeff
    ########################         
    matrixapr         <- matrixint
    
    matrixapr$index <- 1:nrow(matrixapr)          
    mn2       <- merge(matrixapr,MEAN1s, by="sample")
    mn2$M     <- mn2$M - mn2$mean.sample
    matrixapr <- mn2
    matrixapr <- matrixapr[matrixapr$index,]
    nam       <- colnames(matrixint)
    matrixapr <-matrixapr[,nam]
    
    ######################## new  
    #for protein effect and sample for condition effect
    
    matrixeff <- matrixapr
    
    states <- unique(matrixdat$sample)
    for (s in states){
      inds               <- which(matrixapr$sample==states[s])
      matrixeff$M [inds] <- matrixapr$M [inds] -  mean( matrixref$M[which(matrixref$sample==s)]) 
    }
    
    subObjectNorm <- matrixeff

  }   #if (dim(subObjectRef)[1] ==0   
    
   objectNorm    <- rbind(objectNorm,subObjectNorm)
   
} #for (i in )

  return(objectNorm)
}  

  
#----------------------------------------------------------  
objectNorm <-NULL
if         (params$design=="PEP_INTENSITY"){ 
   objectNorm =.effectByProtINT(object,objectRef)
  } else if (params$design=="PEP_RATIO") {     
   objectNorm =.effectByProtRAT(object,objectRef)
}  #if (paramet
#----------------------------------------------------------  
objectNorm = objectNorm%>% filter( !is.na(M))

return(objectNorm)    
   
} ##normEachProt  

####>Revision history<####
# 5.2.4 bug corrected: pAdj.method	fdr and residual.variability are not necessary when normalization.only=yes (IB 30/01/2002)
# 5.2.3 .control function checks parameters for normalisation.method and design (minimal required parameters) (IB 06/01/2020)
# 5.2.2 the R packages are required as the function fread depend on them from now on  (IB 13/11/2019)
# 5.2.1 correction of quanti bug, rename needs dplyr::rename  (IB 06/11/2019)
# 5.2.0 normalization eachProt is able to deel with more than 2 states in both INTENSITY and RATIO designs (IB 17/10/2019)
# 5.1.2 in normalization eachProt, if only State1 remains, protein is not normalized (IB 27/10/2019)  
# 5.1.1 in normalization eachProt changed how to match "states" (IB 10/09/2019)
# 5.1.0 normalization eachProt can be used for several "states" in "PEP_INTENSITY" design (IB 20/08/2019)
# 5.0.7 protein effect in .eachProt normalisation is corrected (IB 09/07/2019)
# 5.0.6 bug correction: protein effect is added instead of subtracted in .normEachProt (IB 08/07/2019)
# 5.0.5 bug correction: parameters instead of params in .normEachProt (IB 03/07/2019)
# 5.0.4 in each Prot normalisation  for INTENSITY design, peptide effect removed from ref (IB 14/06/19)
# 5.0.3 protein effect taken into account in each Prot normalisation  for RATIO design (IB 21/05/19)
# 5.0.2 each Prot normalisation ok for RATIO design (IB 20/05/19)
# 5.0.1 modification how to match proteinId in table.txt and tableRef.tx (IB 17/05/19)
# 5.0.0 normalization eachProt has a new algorithm (quantile without outliers and protein effect) (IB 17/05/19)
# 4.5.3 estimates of sample in .lmProteinPEP_INTENSITY are no longer inverted (IB 15/04/19)  
# 4.5.2 in .analysisDiff the merge of technical replicates is nicer (IB 25/03/19) 
# 4.5.1 if data contains only one protein, function .myCorrectionCI is modified in order to get it into account (IB 15/03/19)
# 4.5.0 if residual.variability=biological and technical replicates existe, these are summarized (mean) and only one pseudo technical replicate is taken into account for the quantification step (.analysisDiff function)(IB 24/01/19)
# 4.4.0 normalization.ref.test has now 4 items, 1 and 2 for normalisation of tableRef and 3 and 4 for correction of data with tableRef (IB 29/11/2018)
# 4.3.1 correct .control if normalization.ref.test doesn't existe (IB 29/11/2018)
# 4.3.0 change Licence from GPL to CeCILL, add dataRef in .control arguments, add in .control function:  Control the correct number and type of files, Control the correct type of the normalization.ref.test parameter (16/11/18)
# 4.2.2 correction of the ic and correction of the normalization in prot (29/06/18)
# 4.2.1 correction of the outlier detection (15/06/18)
# 4.2.0 add the function .normEachProt (12/06/18)
# 4.1.13 add a control function to check the parameters, data and normalization protein (12/06/18)
# 4.1.12 return all the values used in the normalization and correct the bug when there is too much fractions (11/06/18)
# 4.1.11 the previous modification is now also applied in the PEP_RATIO design (31/05/18)
# 4.1.10 proteins observed only on one replicate have now NA as p.values and sd (29/05/18)
# 4.1.9 correct the CI only if the p.value is and return the difference to the mean as the output of the normalization (25/05/18)
# 4.1.8 modification of the model for pep_intensity : Abundance~condition+peptide (07/03/18)
# 4.1.7 correction CI and outlier filtering improved : remove in booth states (08/01/18)
# 4.1.6 Add a progress bar in .analysisDiff and correct a bug in the normalization function (11/12/17)
# 4.1.5 Normalization function are more logical and redable (07/12/17)
# 4.1.4 recover the correct values in resultPep (the outliers) (04/12/17)
# 4.1.3 now the outlier detection is different between the 2 designs PEP_INTENSITY and PEP_RATIO (30/11/17)
# 4.1.2 if the length of the argument of .numOrder is 1 return 1 (27/11/17)
# 4.1.1 new sort function : .numOrder (24/11/17)
# 4.1.0 Add a replicate effect when the design is biological (23/11/17)
# 4.0.9 add the names of the runs in percentageNAsample.txt (10/11/17)
# 4.0.8 Change the names of te designs LABELFREE & LABELED to PEP_INTENSITY and SUPERRATIO to PEP_RATIO (08/11/17)
# 4.0.7 Correction of a bug if no outliers were dected and secure the call of the function rename (+ dplyr::) (03/11/17)
# 4.0.6 when a peptide is missing in one condition (LF and LABELED) we remove it in the others conditions + correction of log2 (30/10/17)
# 4.0.5 Aggregation of the multiple observation in LF & LABELED on the intensity and in SR on the ratios (27/10/17)
# 4.0.4 correction of the function .calculRatio (log the intensity only once) (25/10/17)
# 4.0.3 optimisation of the log2ratio computation & correction of a bug in the quantification of the superratio (24/10/17)
# 4.0.2 Correction of the sign of the ratios in LF and correction of the normalisations (20/10/17)
# 4.0.1 the previous code was badly versioned : the debug lines were not commented (18/10/17)
# 4.0.0 Normalisation changed, models on quantification changed and calcul ratio changed  (13/10/17)
# 3.0.12 calculRatio changed, analysisDiff changed and normalisation changed (13/10/17)
# 3.0.12 conditon -> condition (19/09/17)
# 3.0.11 correction of the script for normalisation with specific proteins and correction of the intensity of the conditions in LABELFREE (15/09/17)
# 3.0.10 the output of the normalisation function is now the data plus the bias, the bias is write in a file on the output section (12/09/17)
# 3.0.9 in labelfree design there is allways a protein effect & the function of normalisation now write the coeff of correction in the file allbias_coef.txt (AS 06/09/17)
# 3.0.8 change the names of the data returned by .analysisDiff (AS 05/09/17)
# 3.0.7 add the primary ratios for LABELFREE (AS 05/09/17)
# 3.0.6 minor correction in .lmProteinLABELFREE (AS 31/08/17)
# 3.0.5 change separation of the QSet in resultsPep.txt (AS 31/08/17)
# 3.0.4 add the functions to filter the outliers & change the names of the parametes (ex Design->design) (AS 18/08/17)
# 3.0.3 optimisation of .computeLog2Ratios  (AS 07/07/17)
# 3.0.2 rebuild the correct information of the biological replicate  (AS 31/05/17)
# 3.0.1 correction of the design super-ratio  (AS 19/05/17)
# 3.0.0 new version of the script, a lot of changes  (AS 09/05/17)
# 2.5.0 the new version of the scipt : a lot of changes (AS 11/05/16)
# 2.5.0 rewrite the normalization function to add globalStandards & correct bug in multi SILAC with replicates (AS 11/05/16)
# 2.4.1 Add the mean value of intensity for each peptide in the table (in /result) log2MeanPep.txt. (AS 11/05/16)
# 2.4.0 Correction of bugs in the new function detectpepout (AS 10/05/16)
# 2.3.9 A new function to detect outliers for all designs (detectpepout), description commented in the function (AS 04/05/16)
# 2.3.8 Correction of bugs in the function ".analysis" (AS 28/04/16)
# 2.3.7 Change the function "detectpepout" in LABELFEE design : MA is now calculated on replicates (instead of peptides). (AS 22/04/16)
# 2.3.6 Change the function .analysis to allways (except outliers) return a ratio and a p-value (some time NA) and add commentary (AS 22/04/16)
# 2.3.5 correct the detection of outliers in LABELFREE design  (AS 18/04/16)
# 2.3.4 correct plots in MissingSample and MatCorPlot  (AS 05/04/16)
# 2.3.3 add commentary and debug the loess normalisation + minors bugs corrected (AS 29/09/15)
# 2.3.2 add handling all protein with one or two proteins (ML 25/02/15)
# 2.3.1 add replicate with multisilac (ML 25/02/15)
# 2.3.0 keep matched and one bias multisilac (ML 23/02/15)
# 2.2.9 optimize Label Free (ML 20/02/15)
# 2.2.8 optimize fast code (ML 17/02/15)
# 2.2.7 optimize fast code (ML 12/02/15)
# 2.2.6 optimize fast code (ML 09/02/15)
# 2.2.5 optimize fast code (ML 05/02/15)
# 2.2.4 rename bias for SUPERSILAC (ML 05/01/15)
# 2.2.3 add data.table package fast function and add test matrix correlation and percentage NA  (ML 18/12/14)
# 2.2.2 impute.method LABELFREE  (ML 17/12/14)
# 2.2.1 optimize graph LABELFREE MULTISILAC (ML 05/12/14)
# 2.2.0 optimize graph SILAC SUPERSILAC (ML 03/12/14)
# 2.1.9 add R correlation in the graph, correction aggregation max fraction, add peptideIDS for imputed value  (ML 27/11/14)
# 2.1.8 add method of aggregation for fraction and method of imputation for missing values, and adjust stat descr  (ML 26/11/14)
# 2.1.7 add  desc Fraction and clustering function  (ML 17/11/14)
# 2.1.6 add stat desc Missing Values and Replicate function  (ML 13/11/14)
# 2.1.5 add PCA function  (ML 12/11/14)
# 2.1.4 optimization code R  (ML 06/11/14)
# 2.1.3 character Protein_ID  (ML 05/11/14)
# 2.1.2 add function for new modelisation  (ML 04/11/14)
# 2.1.1 add TechnReplicate  (ML 30/10/14)
# 2.1.0 optimization code R  (ML 09/10/14)
# 2.0.4 correction all protein with unique peptide  (ML 09/10/14)
# 2.0.3 graph MAplot for multisilac  (ML 08/10/14)
# 2.0.2 correction graph quantile and vsn normalization  (ML 09/09/14)
# 2.0.1 correction MAplot  (ML 08/09/14)
# 2.0.0 add noramlisation  (ML 31/07/14)
# 1.1.0 Optimization SD  (ML 19/06/14)
# 1.0.9 remove Rplots.pdf  (ML 16/06/14)
# 1.0.8 add comments  (ML 12/06/14)
# 1.0.7 add function to adjust ResultsDAprot.txt (ML 11/06/14)
# 1.0.6 remove window pdf and add labelfree descrep(ML 10/06/14)
# 1.0.5 add MAplot Label Free  (ML 03/06/14)
# 1.0.4 add MAplot,vsn normalization, global normalization mean and median (ML 05/06/14)
# 1.0.3 Reverse Ratio and add pooled sd  (ML 03/06/14)
# 1.0.2 add packages  (ML 28/05/14)
# 1.0.1 Refreshment, Remove protein with one peptide, output graph file text protein,add function adjustData and Peptide_IDs  (ML 26/05/14)
# 1.0.0 First version (ML 22/05/14)
