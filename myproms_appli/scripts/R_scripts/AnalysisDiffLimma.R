################################################################################
# AnalysisDiffLimma.R         4.1.2                                            #
# Authors: Matthieu Lhotellier & Alexandre Sta (Institut Curie)                #
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

#### Load functions & libraires ####
source(paste(filepath,"FunctionLimma.R",sep=""))

#### Session info ####
sink(file.path("results","sessionInfo.txt"))
print(sessionInfo())
sink()

#### Load data ####
data<-fread(file.path("data","table.txt"),sep="\t",header=TRUE,data.table=FALSE,integer64 = "numeric")
# n = length(unique(data$Protein_ID)) # DEBUG
# nprot = 500 # DEBUG
# set.seed(25042017)
# data2 = data %>% filter(Protein_ID%in% unique(data$Protein_ID)[sample(1:n,nprot)] ) # DEBUG
# data=data2 # DEBUG
# data = data2 %>% full_join( (data %>% filter(Protein_ID=="mySelected protein")))

data = .reshapeData(data) #%>% .recoverReplicate()
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

#### Control ####

.control(data,parameters,normProtein)

#### Missing sample ####
.missingSample(data,step="Brut_")

#### Compute the log2Ratio  ####
data<-.calculLog2Ratio(data,parameters)
historyData$dataRatio = data

#### Normalisation in prot ####
if( ifelse(!is.null(parameters$normalization.ref),parameters$normalization.ref,"notNone") =="none.none" ){
  dataRef <- fread(file.path("data","tableRef.txt"),sep="\t",header=TRUE,data.table=FALSE,integer64 = "numeric")
  dataRef = .reshapeData(dataRef)
  dataRef = .calculLog2Ratio(dataRef,parameters)
  data = .normEachProt( data , dataRef , parameters )
}

#### Normalisation  ####
tmp = .normalizeData(data,parameters,normProtein)
data <-tmp$data
bias = tmp$bias
rm("tmp")
historyData$dataNorm = data


#### Remove outliers  ####
# historyData$rawData %>% select(sample,experiment,proteinId,proteinName,peptide,peptideId) %>%
#  unique %>% mutate(out="nonOut") %>% spread(sample,out) %>% filter(is.na(State1)|is.na(State2))# DEBUG
# break()
tmp = .myFilterOutlier(data,parameters)
data <-tmp$data
outlier = tmp$outlier
rm("tmp")
historyData$dataWithoutOutlier = data

#### Quantification  ####
data <- .analysisDiff(data,parameters)
historyData$dataQuanti = data

################################################################################
#################################### Output ####################################
################################################################################
################################################################################
#### Bias output ####
if(!is.null(bias))
{
  toWrite = bias %>% (base::t)
  write.table(toWrite,paste("results/allbias_coef.txt"),row.names = TRUE,col.names = FALSE ,sep="\t",quote=FALSE)
}
#### ResultsPep.txt (outliers here) ####
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

#### Build ResultsDAProt.txt ####
# estimate -> log2FC_ for ratios & log2Mean for primary ratios
resultsDAProt = data %>% dplyr::rename(ProteinID=proteinId) %>% .mySpread(value = "estimate", newName = "Log2FC_")
# Change the log2FC into log2Mean for the columns of primary ratios
if(parameters$design=="PEP_INTENSITY"){
  tmp = base::names( resultsDAProt ) %>% 
    base::strsplit( split = "_" ) %>% 
    base::lapply( FUN = function( x ) x[ 2 ] ) %>% 
    ( base::unlist ) 
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

#### Correlation matrix ####
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
  geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                     midpoint = 0, limit = c(-1,1), space = "Lab",
                                     name="Correlation : pairwise.complete.obs") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1)) +
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
  geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                     midpoint = 0, limit = c(-1,1), space = "Lab",
                                     name="Correlation : pairwise.complete.obs") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)
ggsave("results/graph/Afterematrixcorrelation.jpeg",width=10,height=10)

#### Density & boxplot ####
# Differents cases of designs
# if(parameters$design!="LABELFREE"){
#   tmpNonNorm = historyData$dataRatio %>% .MAToIntensity
#   tmpNorm = historyData$dataWithoutOutlier %>% .MAToIntensity
# }else{
  tmpNonNorm = historyData$dataRatio %>% dplyr::rename(intensity=M)
  tmpNorm = historyData$dataWithoutOutlier %>% dplyr::rename(intensity=M)
# }
# Before normalization
tmp = tmpNonNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)

tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE) %>% 
  ggplot(aes(intensity,color=run)) + geom_density() + labs(x = "Log2 intensity")
ggsave("results/graph/Beforealldensity.jpeg",width=10,height=10)
a = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] ) 
a %>% ggplot(aes(run,intensity,color=sample)) + geom_boxplot() + labs(x = "Sample" , y = "Log2 intensity") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("results/graph/Beforeallpeptide.jpeg",width=10,height=10)

# After normalization
tmp = tmpNorm %>% select(proteinId,peptide,experiment,replicate,repTech,sample,intensity)

tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE) %>% 
  ggplot(aes(intensity,color=run)) + geom_density() + labs(x = "Log2 intensity")
ggsave("results/graph/Afteralldensity.jpeg",width=10,height=10)
a = tmp %>% unite(run,experiment,sample,replicate,sep=".",remove=FALSE)
a$run = factor( a$run , levels= unique(a$run)[.numOrder(unique(a$run))] ) 
a %>% ggplot(aes(run,intensity,color=sample)) + geom_boxplot() + labs(x = "Sample" , y = "Log2 intensity") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("results/graph/Afterallpeptide.jpeg",width=10,height=10)

#### Distribution p.values ####
tmp = historyData$dataQuanti %>%  group_by(term) %>% nest
for(i in 1:(dim(tmp)[1]) ){
  title = tmp$term[i]
  path = paste0("results/graph/distribPValue_",title,".jpeg")
  p = tmp$data[[i]] %>% filter(!is.na(p.value)) %>%
    ggplot(aes(p.value)) + geom_histogram(aes(y=..density..),binwidth = 0.01,fill="white",color="black") +
    coord_cartesian(xlim = c(0,1)) + geom_density(color="black") + ggtitle(title)
  ggsave(plot = p,path,width=10,height=10)
}

#### Distribution the fold change ####
tmp = historyData$dataQuanti %>%  group_by(term) %>% nest
for(i in 1:(dim(tmp)[1]) ){
  title = tmp$term[i]
  path = paste0("results/graph/distriblog2FC_",title,".jpeg")
  nProt = tmp$data[[i]] %>% filter(!is.na(estimate)) %>% .$estimate %>% length
  p = tmp$data[[i]] %>% filter(!is.na(estimate)) %>%
    ggplot(aes(estimate)) + geom_histogram(aes(y=..density..),bins = min(nProt,100),fill="white",color="black") +
    geom_density(color="black") + ggtitle(title)
  ggsave(plot = p,path,width=10,height=10)
}

#### Write the design ####
write.table(parameters$design,paste("results/design.txt"),row.names = FALSE,col.names = FALSE,sep="\t",quote=FALSE)

#### End ####

print("End of the quantification")

####>Revision history<####
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
