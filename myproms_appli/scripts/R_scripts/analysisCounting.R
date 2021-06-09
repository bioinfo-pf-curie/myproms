################################################################################
# analysisCounting.R         1.1.1                                             #
# Authors: A. Sta, P. Poullet (Institut Curie)                                 #
# Contact: myproms@curie.fr                                                    #
# Description : Computes SSPA                                                  #
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

# Wokrdirectory
#setwd("")
repertory = "" # Depository repertory, empty or ends by a "/" please

#####################################################
# Library 
#####################################################
library(reshape2)
library(ggplot2)
library(tidyverse)
library(MASS)
library(data.table)

#####################################################
# Load Data
#####################################################
data = read.table("data/table.txt",head=TRUE) # WARNING: only integer values are ollowed

#####################################################
# Change the shape of the data
#####################################################
data$Protein = as.character(data$Protein)
#Data0 = melt(data,id="Protein") # Short to long : one line per count
Data0 = reshape2::melt(data,id="Protein") # Short to long : one line per count
# <- WARNING: reshape2::melt is obsolete: Find alternative (PP 02/07/20)
Data0$variable = as.character(Data0$variable) # factor -> character
Data0$Protein = as.character(Data0$Protein) # factor -> character

Data0$grp = unlist(lapply(Data0$variable,function(b)unlist(strsplit(b,split="\\."))[1])) # Groups information
Data0$sample = unlist(lapply(Data0$variable,function(b)unlist(strsplit(b,split="\\."))[2])) # Sample information

#####################################################
# Detect and exclude outliers
#####################################################
cat("\n######\nRemoving outliers...\n")
remove_outliers <- function(x) {
	qnt <- quantile(x$value, probs=c(.25, .75), na.rm = TRUE)
	H <- 1.5 * IQR(x$value, na.rm = TRUE)
	y <- x
	y$value[x$value < (qnt[1] - H)] <- NA
	y$value[x$value > (qnt[2] + H)] <- NA
	return(y)
}
n = length(unique(Data0$Protein))
c=1
percent = 0 # To print the progression of the analysis
percent.old = -1  # To print the progression of the analysis
DataFrag <- list()
for(i in unique(Data0$Protein)) # Loop on each protein
{
	prot <- Data0[Data0$Protein == i,]
	for(j in unique(prot$grp)) # Loop on each group
	{  
		grp <- prot[prot$grp == j,]
		DataFrag <- c(DataFrag, list( remove_outliers(grp) ))	
	}

  # Progression informations (plot information each %)
  c=c+1
  percent.old = percent
  percent = round(c/n,2)*100
  if(percent.old!=percent)
  {
    cat(paste("######\n",percent," %\n",sep=""))
  }
}
Data <- as.data.frame(rbindlist(DataFrag))
outliers <- Data[is.na(Data$value),c(1,2)]
colnames(outliers)[2] <- "BioRep"
Data <- na.omit(Data)
cat(paste(dim(outliers)[1]," outliers found.\n\n",sep=""))

#####################################################
# Create the output table
#####################################################
data.out = data.frame(Protein = data$Protein) 
data.out$couple = "0"  # Couple selected 
data.out$effect = "-1" # Effect of the second second group compared to the first
data.out$pvalue = "2"  # p-value of the associated efect
data.out[,paste(unique(Data$grp),"mean",sep=".")] = "-1" # mean of each groups
data.out$bestDelta = 0
data.out$model = "a"

#####################################################
# Compute the univariate analysis
#####################################################
cat("\n######\nComputing analysis...\n")
testGlmNb <- function(dt) {
	result <- tryCatch(
		{
			summary(suppressWarnings(glm.nb(as.numeric(value) ~ step,data = dt)))
		},
		error=function(err) {
			return ("error")
		},
		warning = function(war) {
			#return ("***WARNIGNS")
		},
		finally = {
			# do nothing
		}
	)
	return (result)
}

n = length(unique(Data$Protein))
#g = length(unique(Data$grp))
j=1
percent = 0 # To print the progression of the analysis
percent.old = -1  # To print the progression of the analysis
for(i in unique(Data$Protein)) # Loop on each protein for separate analysis
{
#cat(paste(i,"\n",sep="")) # DEBUG <----------------------------------------------------------------
  dt = Data[Data$Protein == i,] # Current data
  dimnames(dt)[[1]] = seq(1:length( dimnames(dt)[[1]] )) # Change the numbering of the data.frame dt
  
  # Compute the mean of each groups, order them and chose the higer steep (ordinate difference)
   mean = unlist(lapply(unique(dt$grp),function(b)mean(dt$value[dt$grp==b]))) # Compute the mean
   names(mean) = unique(dt$grp) # Name the mean with the correct names of the groups
   mean = mean[order(mean,decreasing=TRUE)] # Order the means computed
   df = diff(mean) # Compute the height of the steps
   nn = length(mean) # number of groups
   names(df) = paste(names(mean)[1:(nn-1)],names(df),sep="/") # Name the steps
   cpl.selected = names(which.max(abs(df))) # Chose the best one
  
  # Compute best delta % (PP)
  data.out$bestDelta[data.out$Protein==i] = 100 * max(abs(df)) / mean[1] # negative.binomial & poisson
  #data.out$bestDelta[data.out$Protein==i] = 100 * max(abs(df)) / (mean[1]-mean[nn]) # gaussian
  
  # compute models 
   cp1 = unlist(strsplit(cpl.selected,split="/"))[1] # Extract the first group of the best couple selected
   cp2 = unlist(strsplit(cpl.selected,split="/"))[2] # Extract the second group of the best couple selected
   dt$step = ifelse(dt$grp %in% names(mean)[which(names(mean)==cp2):length(mean)],"right","left" ) # Affect each observation to the corresponding step (only 2 step)
   dt$order =  unlist(lapply(dt$grp,function(b) which(b==names(mean))))
   dt$mean = mean[dt$order]
   # Plot the values if you want to debug the code
    #dt.plot = dt[,c("value","order","grp","step")]
    #extend.dt = data.frame(value=mean,order=1:length(mean),grp=rep("mean",length(mean)),step="left")
    #dt.plot = rbind(dt.plot,extend.dt)
    #p = ggplot(dt.plot,aes(x=order,y=value,color=grp,shape=step)) + geom_point() + ylim(c(0,max(dt$value)))
    #print(p)
    #scan()
   #cpl.selected = paste(c(cp1,cp2)[order(c(cp1,cp2))],collapse="/") # Order the couple selected by alphabetical order like : grp1/grp2
   # Add one observation if a step is empty : it is like adding one sample where the protein have been seen once
			if(sum(as.numeric(dt$value[dt$step=="left"]))==0){
					dt[dim(dt)[1]+1,c("grp","value","step")] = c(cp1,1,"left") 
			}
			if(sum(as.numeric(dt$value[dt$step=="right"]))==0){
					dt[dim(dt)[1]+1,c("grp","value","step")] = c(cp2,1,"right")
			}
   # Compute a Poisson model : model of counting observation, the effect correspond to the effect of the second cluster (alphabetically ordered previously)
   # If the variance of the observations on each conditions are 0 compute a poisson model if not compute a negative binomial model because negative binomial 
   # model can not be computed between conditions with variance=0
			# Update 28/03/20 (PP): Try negative binomial first. If fails fall back to poisson
 			
			p1 = testGlmNb(dt)
			if (class(p1)[1]=="character") { # glm.nb failed => switch to poisson
				p1 = summary(glm(as.numeric(value) ~ step,data = dt,family = "poisson"))
				data.out$model[data.out$Protein==i] = "poisson"
			} else { # success
				data.out$model[data.out$Protein==i] = "negative.binomial"
			}

  # Compute the effect of the difference of groups 
   effect = p1$coefficients[2,"Estimate"] # Effect of grp2 
   pvalue = p1$coefficients[2,"Pr(>|z|)"] # Associated p-value (t-test !=0)
  # Fill the data
   data.out$pvalue[data.out$Protein==i] = pvalue
   left = paste(unique(dt$grp[dt$step=="left"]),collapse=" + ")
   right = paste(paste(unique(dt$grp[dt$step=="right"]),collapse=" + "))
   data.out$couple[data.out$Protein==i] = paste(left,right,sep=" / ")
   data.out$effect[data.out$Protein==i] = effect
   data.out[data.out$Protein==i, paste(names(mean),"mean",sep=".")  ] = mean 
  
  # Progression informations (plot information each %)
   j=j+1
   percent.old = percent
   percent = round(j/n,2)*100
   if(percent.old!=percent)
   {
     cat(paste("######\n",percent," %\n",sep=""))
   }
}

data.out$pvalue = as.numeric(data.out$pvalue) # character -> numeric
data.out$pvalue_adjusted = p.adjust(as.numeric(data.out$pvalue),method="fdr") # compute the fdr correction to the p-value (multiple testing)
data.out$pvalue_adjusted = as.numeric(data.out$pvalue_adjusted) # character -> numeric

#####################################################
# Write the results
#####################################################
file = paste(repertory,"results/outTable.txt",sep="")
#write.table((data.out %>% dplyr::select(-model)),file,sep="\t",row.names = FALSE) # w/o model
write.table(data.out,file,sep="\t",row.names = FALSE,quote=FALSE) # w model
if (dim(outliers)[1] > 0) {
	outlierFile = paste(repertory,"results/outliers.txt",sep="")
	write.table(outliers,outlierFile,sep="\t",row.names = FALSE,quote=FALSE) 
}


#####################################################
# Visualisation (depend on the design of the study -> may change with the analysis)
#####################################################
plots = !TRUE
if(plots == TRUE)
{
  # Reorder the data to fill the object "data"
  data.out = data.out[order(data.out$Protein),]
  data = data[order(data$Protein),]
  # Fill the object "data"
   data$pvalue = data.out$pvalue
   data$pvalue_adjusted = data.out$pvalue_adjusted
   data$couple = data.out$couple
   data$effect = data.out$effect
   data[,paste(unique(Data$grp),"mean",sep=".")] = data.out[,paste(unique(Data$grp),"mean",sep=".")] 
   
  # Plot Volcano plot
   data  %>% ggplot(
     aes(
       as.numeric(effect),
       -log10(as.numeric(pvalue_adjusted)),
       color=couple
       )
     ) + geom_point()
  # Plot boxplot
   data$pvalue = as.numeric(data$pvalue)
   data  = data[order(data$pvalue,decreasing = FALSE),]
   for(j in data$Protein) 
   {
     i = which(data$Protein==j)
     dt = data[i,] %>% dplyr::select(-pvalue,-pvalue_adjusted) %>% melt %>% 
       mutate(variable=as.character(variable)) %>% 
       separate(variable,c("grp","toRemove"),sep="\\.",remove=FALSE) %>% 
       dplyr::select(-toRemove)
     # Color
     cp = data[i,]$couple
     A = cp %>% strsplit(split=" / ") %>% unlist %>% .[1] %>% strsplit(split=" \\+ ") %>% unlist()
     B = cp %>% strsplit(split=" / ") %>% unlist %>% .[2] %>% strsplit(split=" \\+ ") %>% unlist()
     dt$A_against_B = "A"
     dt$A_against_B[dt$grp %in% B] = "B"
     # Boxplot 
       p = ggplot(dt, aes(as.factor(grp), value,fill=A_against_B)) + geom_boxplot() + ggtitle(j)
       print(p)
     #boxplot(value~grp,dt,main = j)
     print(paste(" Name : ",data$Protein[i],sep=""))
     print(paste("   p-value : ",data$pvalue[i],sep=""))
     print(paste("   grp : ", data$couple[i],sep=""))
     print(paste("   effect : ", data$effect[i],sep=""))
     print(paste("   ###",sep=""))
     scan()
   }
   

 # dt.1 = data.frame(value = c(1,0,0,0,0,0, 25,26,27,30,26,32), grp = c("g1","g1","g1","g1","g1","g1","g2" ,"g2","g2","g2","g2","g2"))
#  summary(glm(value ~ grp,data = dt.1,family = "poisson"))
}


#####################################################
#####################################################
####>Revision history<####
# 1.1.1 [FEATURE] Added bestDelta (in %) to output table (PP 02/07/20)
# 1.1.0 [FEATURE] Added outlier filtering based on quantile distribution (PP 28/03/20)
# 1.0.6 if the variance in each conditions is 0 compute a poisson model instead of a negative binomial (AS 01/12/17) 
# 1.0.5 p.value corrected by fdr by default now (AS 31/08/17) 
# 1.0.4 change the model from a lm:poisson to a lm:negative Binomial (AS 10/02/17)
# 1.0.3 rename, call the data in data & result in reslut and remove the bug "+ NA"  (AS 04/08/16)
# 1.0.2 Protein column  (AS 01/08/16)
# 1.0.1 Change the test, now test the difference of left step vs right step  (AS 31/05/16)
# 1.0.0 First version (AS 31/05/16)
