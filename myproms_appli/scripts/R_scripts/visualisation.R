################################################################################
# visualisation.R         2.0.2                                                #
# Authors: Alexandre Sta (Institut Curie)                                      #
# Contact: myproms@curie.fr                                                    #
# Visualisation script                                                         #
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
#  Following variables are needed : 
#            step : differents step of the pipeline
#   visualisation : if FALSE no visualisation are done
#   missingsample : if FALSE no missing sample visualisation 
#      DescIDprot : for boxplot and dataProcessPlots

if(visualisation)
{
  
 if(missingsample)
 {
  #Missing Data
  MissingSample(dataVis,step=data.state)
 }
  
 # Descriptive statistics
 DataDescPlot(dataVis,step=normalisation)
 
}

####>Revision history<####
# 2.0.2 the visualisation is now applied to dataVis, the shape of dataVis must be the same as the shape or areaions (AS 29/01/18) 
# 2.0.1 change the call of MissingSample to speed up the code (AS 24/01/18) 
# 2.0.0 the code is now compatible with MSstats >3 (AS 21/11/17) 
# 1.0.1 Fix the bug log2FC, constraint the name of this column to be log2FC and not logFC (AS 29/11/16) 
# 1.0.0 First version (AS 27/07/16)

