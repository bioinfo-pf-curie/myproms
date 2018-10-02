################################################################################
# AffectParametersSwath.R         2.0.2                                        #
# Authors: Alexandre Sta (Institut Curie)                                      #
# Contact: alexandre.sta@curie.fr                                              #
# Parameters of quantiswath.R                                                  #
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

################################################################################
# Fill the variable with the correct parameters
################################################################################

# Name of the normalisation method (must be one of the option of the function MSstats::dataProcess or it will be FALSE)
normalization=as.character(paramC[which(paramC[,1]=="normalization.method"),2])
normalization = ifelse(length(normalization)==0,"equalizeMedians",normalization) # Default equalizeMedians

# Specification of the matrix of the contrasts of the model of the quati. Must be somethig like : State1/State2;State1/State3;State2/State3...
matrix.contrast=as.character(paramC[which(paramC[,1]=="contrasts.matrix"),2])

# Parameter for removing the missing values between the states : TRUE or FALSE
removeMissingValues=as.character(paramC[which(paramC[,1]=="removeMissingValues"),2])
removeMissingValues = ifelse(length(removeMissingValues)==0,TRUE,removeMissingValues) # Default TRUE

# Number of clusters to use for parallelisation
clusters=as.character(paramC[which(paramC[,1]=="clusters"),2])
clusters = ifelse(length(clusters)==0,1,clusters) # Default 1

# Wich subset of feature to use, all, top3...
featureSubset = as.character(paramC[which(paramC[,1]=="featureSubset"),2])
featureSubset = ifelse(length(featureSubset)==0,"all",featureSubset) # Default all

# If top n selected, specify n
n_top_feature = as.character(paramC[which(paramC[,1]=="n_top_feature"),2])
n_top_feature = ifelse(length(n_top_feature)==0,"3",n_top_feature) # Default 3

# Method to summaryze the frature
summaryMethod = as.character(paramC[which(paramC[,1]=="summaryMethod"),2])
summaryMethod = ifelse(length(summaryMethod)==0,"TMP",summaryMethod) # Default TMP

################################################################################
# Default parameters if it is empty
################################################################################
if(length(normalization)==0)
{
 normalization=FALSE
}

if(length(matrix.contrast)==0){
 vect = as.character(unique(areaions$Condition))
 left = matrix(1:length(vect),nrow=length(vect),ncol=length(vect))
 right = t(left)
 matrix.contrast = paste(paste(vect[left[upper.tri(left)]],vect[right[upper.tri(right)]],sep="/"),collapse=";")
}

if(length(removeMissingValues)==0)
{
 removeMissingValues=TRUE
}

################################################################################
# Normalization specification
################################################################################

protNorm=NULL
if(normalization=="globalStandards")
{
 nameStandards =fread(file.path("data","globalStandards.txt"),data.table=FALSE)
 protNorm = as.character(as.matrix(nameStandards))
}




####>Revision history<####
# 2.0.2 add options in dataProcess and debug the cluster option (AS 19/01/18) 
# 2.0.1 add cluster option (AS 15/01/17) 
# 2.0.0 the code is now compatible with MSstats >3 (AS 21/11/17) 
# 1.0.6 change missingvalues b default to TRUE (AS 27/10/16)
# 1.0.5 list of protein to normalize (globalStandards) is character (AS 13/09/16)
# 1.0.4 move the code of affectDefaultValuesSwath.R in this script (AS 18/08/16)
# 1.0.3 add the case if ou don't normalize (AS 16/08/16)
# 1.0.2 add globalStandards code for normalization (AS 12/08/16)
# 1.0.1 equalizeMedias==constant for MSstats v2 (AS 11/08/16)
# 1.0.0 First version (AS 29/07/16)
