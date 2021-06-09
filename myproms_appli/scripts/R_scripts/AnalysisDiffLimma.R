################################################################################
# AnalysisDiffLimma.R         4.5.5                                            #
# Authors: Matthieu Lhotellier & Alexandre Sta & Isabel Brito (Institut Curie) #
# Contact: myproms@curie.fr                                                    #
# Statiscal Methods for protein quantification by mass spectrometry            #
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

#### 1.   LOAD FUNCTIONS & LIBRAIRIES ############
##################################################

  source(paste(filepath,"FunctionLimma.R",sep=""))

  #### Session info ####
  sink(file.path("results","sessionInfo.txt"))
  print(sessionInfo())
  sink()
  

#### 2.   LOAD DATA                   ############
##################################################

  datal<-fread(file.path("data","table.txt"),sep="\t",header=TRUE,data.table=FALSE,integer64 = "numeric")
  # n = length(unique(data$Protein_ID)) # DEBUG
  # nprot = 500 # DEBUG
  # set.seed(25042017)
  # data2 = data %>% filter(Protein_ID%in% unique(data$Protein_ID)[sample(1:n,nprot)] ) # DEBUG
  # data=data2 # DEBUG
  # data = data2 %>% full_join( (data %>% filter(Protein_ID=="mySelected protein")))
  
  data = .reshapeData(datal) #%>% .recoverReplicate()
  historyData = list(rawData = data)

  #### Load parameters ####
  parameters=fread(file.path("data","param_char.txt"),header=FALSE,sep="\t",data.table=FALSE) %>% spread(V1,V2) # read the parameters
  # parameters$design="PEP_RATIO" # DEBUG

  #### Load normProtein ####
  if(sum(list.files("data/")=="normProtein.txt")){
    normProtein=fread(file.path("data","normProtein.txt"),header=TRUE,sep="\t",data.table=FALSE)
  }else{
    normProtein = NULL
  }

  #### Load dataRef ####
  if( !( is.null(parameters$normalization.ref.test)  ) ){
   dataRefl <- fread(file.path("data","tableRef.txt"),sep="\t",header=TRUE,data.table=FALSE,integer64 = "numeric")
   dataRef = .reshapeData(dataRefl)
   historyDataRef = list(rawDataRef = dataRef) 
  }else{
    dataRef = NULL
  }  

  
#### 3.   PREPROCESS DATA             ############
##################################################
  
  #### Control ####
  
  .control(data,parameters,normProtein,dataRef)
  
  #### Missing sample ####
  .missingSample(data,step="Brut_")
  
  #### Compute the log2Ratio  ####
  data<-.calculLog2Ratio(data,parameters)
  historyData$dataRatio = data
 
  
  if (is.null(parameters$normalization.only)){  
  #------------------------------------------------- 
    
    #### 4.A   TECHNICAL NORMALIZATION     ###########
    ##################################################
    
    #### Normalisation  ####
    tmp = .normalizeData(data,parameters,normProtein)
    data <-tmp$data
    bias = tmp$bias
    rm("tmp")
    historyData$dataNorm = data

    
    #### 5.A   NORMALIZATION BY PROT       ###########
    ##################################################
    
    #### Normalisation in prot ####
    if (!is.null(dataRef)) {
      
      dataRef = .calculLog2Ratio(dataRef,parameters)
      historyDataRef$dataRefRatio = dataRef      
      
      parametersTmp = parameters
      parametersTmp$normalization.method = paste(unlist(strsplit(parameters$normalization.ref.test, split="\\."))[1:2],collapse=".")
      
      tmp =     .normalizeData(dataRef,parametersTmp,normProtein)
      dataRef = tmp$data 
      historyDataRef$dataRefNorm = dataRef       
      
      data = .normEachProt( data , dataRef , parameters )
      historyData$dataNormEachProt = data
      
      if ( prod(is.nan(data$M))==1  | prod(is.na( data$M))==1 ) {
        
        parametersTmp$normalization.method="median.none"
        write.table(parametersTmp$normalization.method,paste("results/new_params.txt"),row.names = FALSE,col.names = FALSE)
        
        dataRef = historyDataRef$dataRefRatio      
        tmp =     .normalizeData(dataRef,parametersTmp,normProtein)
        dataRef = tmp$data 
        historyDataRef$dataRefNorm = dataRef    
        
        tmp = historyData$dataNorm 
        data = tmp
        data = .normEachProt( data , dataRef , parameters )
        historyData$dataNormEachProt = data
      }
    }
    
    #### 6.A   REMOVE OUTLIERS             ###########
    ##################################################
    
  	# historyData$rawData %>% select(sample,experiment,proteinId,proteinName,peptide,peptideId) %>%
  	#  unique %>% mutate(out="nonOut") %>% spread(sample,out) %>% filter(is.na(State1)|is.na(State2))# DEBUG
  	# break()
  	tmp = .myFilterOutlier(data,parameters)
  	data <-tmp$data
  	outlier = tmp$outlier
  	rm("tmp")
  	historyData$dataWithoutOutlier = data

  	
  	#### 7.A   QUANTIFICATION              ############
  	##################################################
  
  	if (nrow(   datal)!=sum(   datal$Protein_Validity)) {
  	  tmp1 = datal %>% dplyr::mutate(proteinId =Protein_ID, 
  	                                 peptide   =Peptide, 
  	                                 experiment=Experiment,
  	                                 replicate =Replicate,
  	                                 repTech   =Technical_Replicate,
  	                                 sample    =Sample,
  	                                 quantifSet=Quantif_Set, 
  	                                 peptideId =as.character(Peptide_IDs),
  	                                 intensity =Value )
  	  tmp2 =  tmp1 %>% select(peptideId,validity=Protein_Validity)
  	  tmp  =  inner_join(historyData$dataWithoutOutlier, tmp2, by = "peptideId")  %>% 
  	     mutate(validity=factor(validity), intensity=M) %>% 
  	     filter(validity==1)
  	  
  	  if (nrow(tmp)!=0) {
  	  data <- .analysisDiff(tmp,parameters)
  	  historyData$dataQuanti = data 
  	  } else {print("All Protein Validity are zeros")}
  	} else {
  	#### Quantification  ####
  	data <- .analysisDiff(data,parameters)
  	historyData$dataQuanti = data  
  	}

    } else {
	#-------------------------------------------------
  
  	#### 4.B   TECHNICAL NORMALIZATION     ############
  	##################################################
  
  	#### Normalisation  ####
  	tmp = .normalizeData(data,parameters,normProtein)
  	data <-tmp$data
 	  bias = tmp$bias
  	rm("tmp")
  	historyData$dataNorm = data
  
  	#### 6.B   REMOVE OUTLIERS             ###########
  	##################################################
  	
 	# historyData$rawData %>% select(sample,experiment,proteinId,proteinName,peptide,peptideId) %>%
 	#  unique %>% mutate(out="nonOut") %>% spread(sample,out) %>% filter(is.na(State1)|is.na(State2))# DEBUG
 	# break()
  	tmp = .myFilterOutlier(data,parameters)
 	  data <-tmp$data
 	  outlier = tmp$outlier
 	  rm("tmp")
 	  historyData$dataWithoutOutlier = data
  
	}
	#-------------------------------------------------
  
 
#### 8.   OUTPUT                      ############
##################################################

  ################################################################################
  #################################### Output ####################################
  ################################################################################
  ################################################################################

#### 8.1.   #### Bias output ####
################
  if(!is.null(bias)){
    toWrite = bias %>% (base::t)
    write.table(toWrite,paste("results/allbias_coef.txt"),row.names = TRUE,col.names = FALSE ,sep="\t",quote=FALSE)
  }


#### 8.2.  #### ResultsPep.txt (outliers here) ####
###############  
  if (is.null(parameters$normalization.only)){
  #------------------------------------------------- 
  	if ( !is.null(dataRef) ){
  		resultsPep = historyData$dataNormEachProt %>% dplyr::full_join(outlier,by=names(.)) %>% 
    	dplyr::rename( Condition = sample,
                 ProteinID = proteinId, 
                 Peptide = peptide,
                 Experiment = experiment,
                 log2Measure = M,
                 PeptideId = peptideId
   	   ) %>% dplyr::select(-A,-quantifSet)
    } else if ( is.null(dataRef) ){
    	resultsPep = historyData$dataNorm %>% dplyr::full_join(outlier,by=names(.)) %>% 
      	dplyr::rename( Condition = sample,
                   ProteinID = proteinId, 
                   Peptide = peptide,
                   Experiment = experiment,
                   log2Measure = M,
                   PeptideId = peptideId
        ) %>% dplyr::select(-A,-quantifSet)
    }
  order = c("Experiment","Condition","replicate","repTech","ProteinID","normProtein","Peptide","PeptideId","log2Measure","out")
  resultsPep = resultsPep[,order]
  write.table(resultsPep,paste("results/resultsPep.txt"),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)

  } else {                                                                 
  #-------------------------------------------------
    resultsPep = historyData$dataNorm %>% dplyr::full_join(outlier,by=names(.)) %>% 
      dplyr::rename( Condition = sample,
                   ProteinID = proteinId, 
                   Peptide = peptide,
                   Experiment = experiment,
                   log2Measure = M,
                   PeptideId = peptideId
     ) %>% dplyr::select(-A,-quantifSet)
     order = c("Experiment","Condition","replicate","repTech","ProteinID","normProtein","Peptide","PeptideId","log2Measure","out")
     resultsPep = resultsPep[,order]
     write.table(resultsPep,paste("results/resultsPep.txt"),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
  }


#### 8.3. #### Build ResultsDAProt.txt ####
###############  
  if (is.null(parameters$normalization.only)){ 
  #------------------------------------------------- 

  	# estimate -> log2FC_ for ratios & log2Mean for primary ratios
  	resultsDAProt = data %>% dplyr::rename(ProteinID=proteinId) %>% .mySpread(value = "estimate", newName = "Log2FC_")
  	# Change the log2FC into log2Mean for the columns of primary ratios
  	if(parameters$design=="PEP_INTENSITY"){
   	 tmp = base::names( resultsDAProt ) %>% base::strsplit( split = "_" ) %>% base::lapply( FUN = function( x ) x[ 2 ] ) %>%  ( base::unlist ) 
     names( resultsDAProt )[ tmp %in% unique(historyData$dataWithoutOutlier$sample) ]  = tmp %in% unique(historyData$dataWithoutOutlier$sample) %>% tmp[ . ] %>% 
     paste("Log2Mean_",.,sep="")
    }

    # p.value -> pval_
    tmp = data %>% dplyr::rename(ProteinID=proteinId) %>% .mySpread(value="p.value",newName="pVal",sep="_")
    resultsDAProt = resultsDAProt %>% join(.,tmp,"ProteinID")
    # std.error -> SD_
    tmp = data %>% dplyr::rename(ProteinID=proteinId) %>% .mySpread(value="std.error",newName="SD",sep="_")
    resultsDAProt = resultsDAProt %>% join(.,tmp,"ProteinID")
    # ci2.5 -> CI2.5_
    tmp = data %>% dplyr::rename(ProteinID=proteinId) %>% .mySpread(value= "ci2.5",newName= "CI2.5",sep="_")
    resultsDAProt = resultsDAProt %>% join(.,tmp,"ProteinID")
    # ci97.5 -> CI97.5_
    tmp = data %>% dplyr::rename(ProteinID=proteinId) %>% .mySpread(value= "ci97.5",newName= "CI97.5",sep="_")
    resultsDAProt = resultsDAProt %>% join(.,tmp,"ProteinID")
    # Write the table in the file results
    write.table(resultsDAProt,paste("results/ResultsDAProt.txt"),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
  }

  
#### 8.4. #### Correlation matrix ####
###############  
  # Differents cases of designs
  if(parameters$design!="PEP_INTENSITY"){
    tmpNonNorm = historyData$dataRatio %>% .MAToIntensity
    tmpNorm = historyData$dataWithoutOutlier %>% .MAToIntensity
  }else{
    tmpNonNorm = historyData$dataRatio %>% dplyr::rename(intensity=M)
    tmpNorm = historyData$dataWithoutOutlier %>% dplyr::rename(intensity=M)
  }

  # Data before normalization
  # Reshape the data
  tmp = tmpNonNorm %>% dplyr::select(sample,experiment,replicate,repTech,proteinId,peptide,intensity) %>%
    tidyr::unite(run,sample,experiment,replicate,repTech) %>% tidyr::spread(key=run , value=intensity) %>%
    dplyr::select(-proteinId,-peptide) %>% (base::as.matrix)
  # Matrice de correlation
  matrix.correlation = cor(tmp,use="pairwise.complete.obs")
  write.table(matrix.correlation,paste("results/Beforematrixcorrelation.txt"),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)

  # Plot the correlation matrix
  ggplot(data = melt(matrix.correlation), aes(x=Var1, y=Var2, fill=value)) +
         geom_tile() + 
         scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                              midpoint = 0, limit = c(-1,1), space = "Lab",
                              name="Correlation : pairwise.complete.obs") +
         theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
         coord_fixed() +
         geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)
  ggsave("results/graph/Beforematrixcorrelation.jpeg",width=10,height=10)

  # Data After normalization
  # Reshape the data
  tmp = tmpNorm %>% dplyr::select(sample,experiment,replicate,repTech,proteinId,peptide,intensity) %>%
    tidyr::unite(run,sample,experiment,replicate,repTech) %>% tidyr::spread(key=run , value=intensity) %>%
    dplyr::select(-proteinId,-peptide) %>% (base::as.matrix)
  # Matrice de correlation
  matrix.correlation = cor(tmp,use="pairwise.complete.obs")
  write.table(matrix.correlation,paste("results/Aftermatrixcorrelation.txt"),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)

  # Plot the correlation matrix
  ggplot(data = melt(matrix.correlation), aes(x=Var1, y=Var2, fill=value)) +
         geom_tile() + 
         scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                              midpoint = 0, limit = c(-1,1), space = "Lab",
                              name="Correlation : pairwise.complete.obs") +
         theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
         coord_fixed() +
         geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)
  ggsave("results/graph/Afterematrixcorrelation.jpeg",width=10,height=10)


#### 8.5. #### Density & boxplot ####
###############  
  # Differents cases of designs
  # if(parameters$design!="LABELFREE"){
  #   tmpNonNorm = historyData$dataRatio %>% .MAToIntensity
  #   tmpNorm = historyData$dataWithoutOutlier %>% .MAToIntensity
  # }else{ }

  #################################################################################
  divisors <- function(x){
    #  Vector of numberes to test against
    y <- seq_len(x)
    #  Modulo division. If remainder is 0 that number is a divisor of x so return it
    y[ x%%y == 0 ]
  }  
  
  boxplots = function(a,params,graphname){
    locate <- function(x, targets) {
    results <- lapply(targets, function(target) which(x == target))
    results
    }
    nbbxfixe=30
    if( !( is.null(parameters$nb.boxplots)  ) ){ 
      nbbx=as.numeric(parameters$nb.boxplots)
    }else{ nbbx= nbbxfixe }
    print(nbbx)
    pnbgp=divisors(dim(table(a$sample))*dim(table(a$experiment)))
    bnbgp=dim(table(a$run))/pnbgp
    mnbgp=pnbgp[which(abs(bnbgp-(nbbx+1)) ==min(abs(bnbgp-(nbbx+1))))]
    
    aux1 =  unlist(   strsplit(names(table(a$run)),".r"))           
    aux2 =  aux1 [seq(1, length( aux1    )  ,2    )]                
    aux3 =  unique(aux2)                                             
    aux4 =  split(aux3, rep(c(1:mnbgp),each=length(aux3)/mnbgp) )   
    aux5 =  paste(a$experiment,a$sample,sep=".")                     
    aux6 =  lapply(aux4, function(x){  unlist(locate(aux5,x)) } )
    
    for (i in 1:mnbgp) { print(i)  
      ax = a[aux6[[i]],]
      ax %>% 
        ggplot(aes(run,intensity,color=sample)) + 
        geom_boxplot() + 
        labs(x = "Sample" , y = "Log2 intensity") + 
        scale_y_continuous(limits=range(a$intensity)) +
        theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
      if (mnbgp==1) {
        ggsave(paste0(graphname,".jpeg"),width=10,height=10)
      }else{
      ggsave(paste0(graphname,i,".jpeg"),width=10,height=10)
      }  
    }  # for (i in 
  } 
  
  boxplots_validity = function(a,params,graphname){
    locate <- function(x, targets) {
    results <- lapply(targets, function(target) which(x == target))
    results
    }
    nbbxfixe=30
    if( !( is.null(parameters$nb.boxplots)  ) ){ 
      nbbx=as.numeric(parameters$nb.boxplots)
    }else{ nbbx= nbbxfixe }
    print(nbbx)
    pnbgp=divisors(dim(table(a$sample))*dim(table(a$experiment)))
    bnbgp=dim(table(a$run))/pnbgp
    mnbgp=pnbgp[which(abs(bnbgp-(nbbx+1)) ==min(abs(bnbgp-(nbbx+1))))]
    
    aux1 =  unlist(   strsplit(names(table(a$run)),".r"))  
    aux2 =  aux1 [seq(1, length( aux1    )  ,2    )]
    aux3 =  unique(aux2)
    aux4 =  split(aux3, rep(c(1:mnbgp),each=length(aux3)/mnbgp) )  
    aux5 =  paste(a$experiment,a$sample,sep=".")
    aux6 =  lapply(aux4, function(x){  unlist(locate(aux5,x)) } )
    
    for (i in 1:mnbgp) { print(i)  
      ax = a[aux6[[i]],]
      ax %>% 
        ggplot(aes(run,intensity,color=sample, fill = validity)) + 
        geom_boxplot() + 
        scale_fill_brewer(palette="Dark2") +
        labs(x = "Sample" , y = "Log2 intensity") + 
        scale_y_continuous(limits=range(a$intensity)) +
        theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
      ggsave(paste0(graphname,i,".jpeg"),width=10,height=10)
    }  # for (i in 
  } 
  #################################################################################  
  
    tmpNonNorm    = historyData$dataRatio          %>% dplyr::rename(intensity=M)
    tmpNorm       = historyData$dataWithoutOutlier %>% dplyr::rename(intensity=M)
 
    tmpDATNonNorm = historyData$dataRatio          %>% dplyr::rename(intensity=M) 
    tmpDATNorm    = historyData$dataNorm           %>% dplyr::rename(intensity=M) 
    tmpDATfin     = historyData$dataWithoutOutlier %>% dplyr::rename(intensity=M)
   
  # Before normalization
   tmp = tmpNonNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)
   tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE) %>% 
   ggplot(aes(intensity,color=run)) + 
          geom_density() + 
          labs(x = "Log2 intensity")
   ggsave("results/graph/Beforealldensity.jpeg",width=10,height=10)

   tmp = tmpDATNonNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)
   a = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
   a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] ) 
   
   graphname = "results/graph/Beforeallpeptide"
   boxplots(a,params,graphname)
   #-------------------------------------------
   #a %>% 
   #ggplot(aes(run,intensity,color=sample)) + 
   #       geom_boxplot() + 
   #       labs(x = "Sample" , y = "Log2 intensity") + 
   #       theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
   #ggsave("results/graph/Beforeallpeptide.jpeg",width=10,height=10)
   
  # After normalization
  tmp = tmpNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)
  tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE) %>% 
  ggplot(aes(intensity,color=run)) + 
         geom_density() + 
         labs(x = "Log2 intensity")
  ggsave("results/graph/Afteralldensity.jpeg",width=10,height=10)

  tmp = tmpDATNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)
  a = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
  a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] ) 
  
  graphname = "results/graph/Afternormallpeptide"
  boxplots(a,params,graphname)
  #-------------------------------------------
  #a %>% 
  #ggplot(aes(run,intensity,color=sample)) + 
  #       geom_boxplot() + 
  #       labs(x = "Sample" , y = "Log2 intensity") + 
  #       theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
  #ggsave("results/graph/Afternormallpeptide.jpeg",width=10,height=10)

  tmp = tmpDATfin %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)
  a = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
  a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] ) 
  
  graphname = "results/graph/Afterallpeptide"
  boxplots(a,params,graphname)
  #-------------------------------------------
  #a %>% 
  #ggplot(aes(run,intensity,color=sample)) + 
  #       geom_boxplot() + 
  #       labs(x = "Sample" , y = "Log2 intensity") + 
  #       theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
  #ggsave("results/graph/Afterallpeptide.jpeg",width=10,height=10)

  
  if (is.null(parameters$normalization.only)){ 
  #------------------------------------------------- 
    if (!is.null(dataRef)) {
     tmpREFNonNorm = historyDataRef$dataRefRatio %>% dplyr::rename(intensity=M)
     tmpREFNorm    = historyDataRef$dataRefNorm %>% dplyr::rename(intensity=M)
   
     tmp = tmpREFNonNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)
     a = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
     a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] ) 
     
     graphname = "results/graph/BeforeREFallpeptide"
     boxplots(a,params,graphname)
     #-------------------------------------------
     #a %>% 
     #ggplot(aes(run,intensity,color=sample)) + 
     #        geom_boxplot() + 
     #        labs(x = "Sample" , y = "Log2 intensity") + 
     #       theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
     #ggsave("results/graph/BeforeREFallpeptide.jpeg",width=10,height=10)
    
     tmp = tmpREFNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)
     a = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
     a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] ) 
     
     graphname = "results/graph/AfternormREFallpeptide"
     boxplots(a,params,graphname)
     #-------------------------------------------
     #a %>% 
     #ggplot(aes(run,intensity,color=sample)) + 
     #        geom_boxplot() + 
     #        labs(x = "Sample" , y = "Log2 intensity") + 
     #        theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
     #ggsave("results/graph/AfternormREFallpeptide.jpeg",width=10,height=10)
    }
  } 
  
  if (nrow(   datal)!=sum(   datal$Protein_Validity)) {
  #------------------------------------------------- 
    tmp1 = datal %>% dplyr::mutate(proteinId =Protein_ID, 
                                   peptide   =Peptide, 
                                   experiment=Experiment,
                                   replicate =Replicate,
                                   repTech   =Technical_Replicate,
                                   sample    =Sample,
                                   quantifSet=Quantif_Set, 
                                   peptideId =as.character(Peptide_IDs),
                                   intensity =Value )
    tmp2 =  tmp1 %>% select(peptideId,validity=Protein_Validity)
    tmp  =  inner_join(historyData$dataRatio, tmp2, by = "peptideId")  %>% mutate(validity=factor(validity), intensity=M)
    a    = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
    a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] )  
    
    graphname = "results/graph/Beforeallpeptide_validity"
    boxplots_validity(a,params,graphname)
    #-------------------------------------------
    #a %>% 
    #  ggplot(aes(run,intensity,color=sample, fill = validity )) + 
    #  geom_boxplot() + 
    #  scale_fill_brewer(palette="Dark2") +
    #  labs(x = "Sample" , y = "Log2 intensity") + 
    #  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) 
    #ggsave("results/graph/Beforeallpeptide_validity.jpeg",width=10,height=10)
    
    
    tmp  =  inner_join(historyData$dataWithoutOutlier, tmp2, by = "peptideId")  %>% mutate(validity=factor(validity), intensity=M)
    a    = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
    a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] )  
    
    graphname = "results/graph/Afterallpeptide_validity"
    boxplots_validity(a,params,graphname)
    #-------------------------------------------
    #a %>% 
    #  ggplot(aes(run,intensity,color=sample, fill = validity )) + 
    #  geom_boxplot() + 
    # scale_fill_brewer(palette="Dark2") +
    #  labs(x = "Sample" , y = "Log2 intensity") + 
    #  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) 
    #ggsave("results/graph/Afterallpeptide_validity.jpeg",width=10,height=10)
    
    
    tmp  =  inner_join(historyData$dataNorm, tmp2, by = "peptideId")  %>% mutate(validity=factor(validity), intensity=M)
    a    = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
    a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] )  
    
    graphname = "results/graph/Afternormallpeptide_validity"
    boxplots_validity(a,params,graphname)
    #-------------------------------------------
    #a %>% 
    #  ggplot(aes(run,intensity,color=sample, fill = validity )) + 
    #  geom_boxplot() + 
    #  scale_fill_brewer(palette="Dark2") +
    #  labs(x = "Sample" , y = "Log2 intensity") + 
    # theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) 
    #ggsave("results/graph/Afternormallpeptide_validity.jpeg",width=10,height=10)
  }
  

#### 8.6. #### Distribution p.values ####
###############
  if (is.null(parameters$normalization.only)){ 
  #------------------------------------------------- 
  tmp = historyData$dataQuanti %>%  group_by(term) %>% nest
  for(i in 1:(dim(tmp)[1]) ){
    title = tmp$term[i]
    path = paste0("results/graph/distribPValue_",title,".jpeg")
    p = tmp$data[[i]] %>% filter(!is.na(p.value)) %>%
    ggplot(aes(p.value)) + 
           geom_histogram(aes(y=..density..),binwidth = 0.01,fill="white",color="black") +
           coord_cartesian(xlim = c(0,1)) + 
           geom_density(color="black") + 
           ggtitle(title)
    ggsave(plot = p,path,width=10,height=10)
   }
  }


#### 8.7. #### Distribution the fold change ####
###############
  if (is.null(parameters$normalization.only)){ 
  #------------------------------------------------- 
  tmp = historyData$dataQuanti %>%  group_by(term) %>% nest
  for(i in 1:(dim(tmp)[1]) ){
    title  = tmp$term[i]
    path   = paste0("results/graph/distriblog2FC_",title,".jpeg")
    nProt  = tmp$data[[i]] %>% filter(!is.na(estimate)) %>% .$estimate %>% length
    for_outliers_step1 = tmp$data[[i]] %>% filter(!is.na(estimate)) %>%
      dplyr::summarise(outq1 = quantile(estimate, probs=0.25) - 1.5*IQR(estimate),
                     outq3 = quantile(estimate, probs=0.75) + 1.5*IQR(estimate) )
    for_outliers_step2 = tmp$data[[i]] %>% filter(!is.na(estimate)) %>% 
      filter (estimate < for_outliers_step1$outq1 | estimate > for_outliers_step1$outq3 )  %>% .$estimate
    noutliers = for_outliers_step2  %>% length
    routliers = for_outliers_step2  %>% range 
    title = paste(title, paste0("(There are ",noutliers," outliers, min is ", 
                              round(routliers[1],0), ", max is ", round(routliers[2],0),")"), sep= "     ")
  
    p = tmp$data[[i]] %>% filter(!is.na(estimate)) %>% 
     filter (estimate > for_outliers_step1$outq1 & estimate < for_outliers_step1$outq3 ) %>% 
    ggplot(aes(estimate)) + 
           geom_histogram(aes(y=..density..),bins = min(nProt,100),fill="white",color="black") +
           geom_density(color="black") + 
           ggtitle(title)
    ggsave(plot = p,path,width=10,height=10)
   }
  }

#### 8. #### Write the design ####
###############
  write.table(parameters$design,paste("results/design.txt"),row.names = FALSE,col.names = FALSE,sep="\t",quote=FALSE)


#### 9.   END                         ############
##################################################
#### End #### 
#print("End of the quantification")
print("End of the analysis")

####>Revision history<####
# 4.5.5 the nb of graphics is corrected, when the nb of states is high and several graphics are displayed,  (IB 14/09/2020)  
# 4.5.4 bugs on boxplots of peptides when several are displayed and error incompatibility of peptides columns type are corrected (IB 25/08/2020)   
# 4.5.3 If the nb of states is high, several boxplots are displayed (IB 01/05/2020)  
# 4.5.2 Distribution the the fold change is not available, when normalization.only=yes and normalization.method becomes "median.none" in new_params.txt file istead of "none.none"  (IB 30/01/2020) 
# 4.5.1 If Protein_Validity equal to zero, then new graphs Beforeallpeptide_validity.jpg, Afternormallpeptide_validity.jpg, Afterallpeptide_validity.jpg and removed for quantification; distriblog2FC_....jpg allowed when normalization.only is yes (IB 06/01/2020)   
# 4.5.0 this script is able to read a new parameter in the file "param_char.txt" called normalization.only; if it existes (normalization.only	yes), then only the technical normalization is performed (IB 31/12/2019) 
# 4.4.2 if the (technical) normalisation of dataref fails, parameters are changent to none.none and a file named "new_parameters.txt" is created (IB 17/10/2019)
# 4.4.1 rename boxplot graphs: Beforeallpeptide, BeforeREFallpeptide,  Afternormallpeptide,  AfternormREFallpeptide and Afterallpeptide (IB 20/08/2019)
# 4.4.0 rename graphs Afterallpeptide as AftereachProtnoOutDATallpeptide and Beforeallpeptide as BeforeDATallpeptide; create graphs AfternormProtDATallpeptide,  AfternormProtRefallpeptide and BeforeRefallpeptide (IB 20/08/2019)
# 4.3.2 if tableRef doesn't exist, resultsPep uses the dataNorm (IB 05/04/19)
# 4.3.1 resultsPep uses the dataNormEachProt instead of dataNorm (IB 25/03/19)
# 4.3.0 normalization.ref.test has now 4 items, 1 and 2 for normalisation of tableRef and 3 and 4 for correction of data with tableRef (IB 29/11/2018)
# 4.2.0 change Licence from GPL to CeCILL, change order of table.ref reading, create dataRef object (NULL when there isn't table.ref) and add dataRef in .control arguments, remove outliers from data in the graph "distriblog2FC_....jpeg" (16/11/18)
# 4.1.3 invert the order of normalization and normalization by ref and correct the test on the parameter of normalization of ref (28/06/18)
# 4.1.2 correction of the test on normalization.ref (15/06/18)
# 4.1.1 do not normalize in each protein if normalization.ref is none.none (12/06/18)
# 4.1.0 add the normalization in the protein (12/06/18)
# 4.0.13 add a control section to check the parameters, data and normalization protein (12/06/18)
# 4.0.12 minor correction of commetary and the call of functions (08/01/18)
# 4.0.11 add a progress bar to analysisDiff (11/12/17)
# 4.0.10 force the value of intensity to be read as a numerical value instead of a integer (08/12/17)
# 4.0.9 the outlier detection now needs the information of the parameters (30/11/17)
# 4.0.8 correction of the distribution of the FC if there is less than 100 proteins (27/11/17)
# 4.0.7 arrange the graph Afterallpeptide and Beforeallpeptide (24/11/17)
# 4.0.6 Add te distribution of the FC (23/11/17)
# 4.0.5 design becomes PEP_INTENSITY or PEP_RATIO (10/11/17)
# 4.0.4 Reorder resultsPep (27/10/17)
# 4.0.3 Add the design file (27/10/17)
# 4.0.2 Correction of the name of the columns of ResultDAProt for LF and SR (24/10/17)
# 4.0.1 Correction of the sign of the ratios in LF and correction of the normalisations (20/10/17)
# 4.0.0 Normalisation changed, models on quantification changed and calcul ratio changed  (13/10/17)
# 3.0.13 normProtein is a data.frame (19/09/17)
# 3.0.12 load the normProtein with the header (15/09/17)
# 3.0.11 uniformalisation of the names of the runs : condition_experiment_repbio_reptech & write the bias in the output section (12/09/17)
# 3.0.10 split the data in ResultsPep.txt (12/09/17)
# 3.0.9 don't change the names of the column of resultsDAProt.txt in the LABELED design (06/09/17)
# 3.0.8 move the writing of allbias_coef.txt into the function of normalization (06/09/17)
# 3.0.7 change the name of the column of the primary ratios by log2Mean_ in resultsDAProt.txt (AS 05/09/17)
# 3.0.6 don't create the file ResultsProt.txt (AS 31/08/17)
# 3.0.5 remove the "" in the .txt files & change separation of the QSet in resultsPep.txt (AS 31/08/17)
# 3.0.4 correction of newName in .mySpread (AS 30/08/17)
# 3.0.3 remove param numeric file (AS 21/08/17)
# 3.0.2 add the detection and filtering of the outliers & change the names of the parametes (ex Design->design) (AS 18/08/17)
# 3.0.1 rebuild the correct information of the biological replicate  (AS 31/05/17)
# 3.0.0 new version of the script, a lot of changes  (AS 09/05/17)
# 2.2.5 correction of the default value of the bias.correction  (AS 10/11/16)
# 2.2.4 correction of the FDR  (AS 28/04/16)
# 2.2.3 modification of the output "ResultsDAProt.txt" : case of missing value of protein intensity in an experiment (AS 26/04/16)
# 2.2.2 same as 2.3.1PP (AS 05/10/15)
# 2.3.1PP HeatmapPvalue call commented (PP 12/06/15)
# 2.2.1 add graph after normalization and outliers  (ML 23/02/14)
# 2.2.0 add superratio one peptide  (ML 17/02/14)
# 2.1.9 add data.table package fast function  (ML 17/12/14)
# 2.1.8 optimize process graph  (ML 03/12/14)
# 2.1.7 add stat desc  (ML 17/11/14)
# 2.1.6 optimization code R  (ML 05/11/14)
# 2.1.5 correction LabelFree  (ML 04/11/14)
# 2.1.4 add column Experiment  (ML 03/11/14)
# 2.1.3 correction Experiment  (ML 03/11/14)
# 2.1.2 add AllRatioProt  (ML 03/11/14)
# 2.1.1 add TechnReplicate  (ML 30/10/14)
# 2.1.0 optimization code R  (ML 09/10/14)
# 2.0.8 correction all protein with unique peptide  (ML 09/10/14)
# 2.0.7 GraphMAplot for multiSILAC (ML 08/10/14)
# 2.0.6 Correction Peptide_IDs for pipeline (ML 02/10/14)
# 2.0.5 Peptide_IDs for multiSILAC (ML 30/09/14)
# 2.0.4 filepath function and correction multiplex (ML 26/09/14)
# 2.0.3 add analysis multiplex labelled (ML 10/09/14)
# 2.0.2 delete graph afterfirstnorm and files and optimization choice normalization (ML 09/09/14)
# 2.0.1 Optimization time for grubb test and correction LabelFree (ML 08/09/14)
# 2.0.0 Add Normalization (ML 01/08/14)
# 1.1.7 Modification Normalization (ML 24/06/14)
# 1.1.6 Modification Filter peptide matched for ratio (ML 23/06/14)
# 1.1.5 venn diagram (ML 20/06/14)
# 1.1.4 Correction ratio and venn diagram (ML 20/06/14)
# 1.1.3 optimization Peptide_IDs and SD and correction graph file (ML 19/06/14)
# 1.1.2 add the choice to find or not Peptide_IDs and SD (ML 16/06/14)
# 1.1.1 add comments and correction when there is one experiment (ML 12/06/14)
# 1.0.9 adjust ResultsDAplot.txt  (ML 11/06/14)
# 1.0.8 All Output Data log 2  (ML 10/06/14)
# 1.0.7 specify MAplot with Design (ML 06/06/14)
# 1.0.6 add MA plot and argument normalization.method  (ML 05/06/14)
# 1.0.5 Reverse Ratio, Rename Graph File and add sd, pooled sd  (ML 03/06/14)
# 1.0.4 Separate graphics ratio proteins  (ML 28/05/14)
# 1.0.3 Correction values and add Peptide_IDs version 2  (ML 27/05/14)
# 1.0.2 List of protein with one peptide, Refreshment, adjust Data and add Peptide_IDs  (ML 26/05/14)
# 1.0.1 Remove filter same proteins (ML 23/05/14)
# 1.0.0 First version (ML 22/05/14)    
