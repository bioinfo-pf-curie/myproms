################################################################################
# motifEnrichment.R       1.0.0                                                        #
# Authors: Stephane Liva    #
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

rm(list=ls())
require(ggplot2)
require(ggseqlogo)
source("/bioinfo/pipelines/myproms/dev/scripts/R_scripts/sliva/motifx-functions.R")

args<-commandArgs(trailingOnly = TRUE)

fgSeq <- as.character(args[1])
centralRes<-as.character(args[2])
occurence<-as.numeric(args[3])
pvalCutoff<-as.numeric(args[4])
bgValue<-as.character(args[5])
width<-as.numeric(args[6])
seqSize<-(width-1)/2
backNbSeq<-as.numeric(args[7])

if (bgValue == 'random') {
	freqAA<-read.table("randomCompAA.txt", sep="\t",row.names=1)
	colnames(freqAA)<-c("FREQ")
	bgSeq<-replicate(backNbSeq,
		  paste0(c(sample(as.vector(rownames(freqAA)),
			   size=seqSize,
			   prob=as.numeric(freqAA$FREQ),
			   replace=TRUE),
			   sample(centralRes),
			   sample(as.vector(rownames(freqAA)),
			   size=seqSize,
			   prob=as.numeric(freqAA$FREQ),
			   replace=TRUE )),
			   collapse=""))
	write.table(bgSeq,"backgroundRandom.txt",col.names=F,quote=F, row.names=F)
	bg.Seq<-readLines("backgroundRandom.txt")
}
if (bgValue == 'quanti') {
	##BACKGROUND QUANTI/PROJET
	bg.Seq<-readLines("backgroundQuanti.txt")
}
fg.Seq<-readLines(fgSeq)

motifTable<-motifx(fg.Seq, bg.Seq, central.res = centralRes, min.seqs = occurence, pval.cutoff = pvalCutoff)
write.table(motifTable, "motifResult.txt", col.names=T, row.names=T, sep="\t")

####>Revision history<####
# 1.0.0 first version (SL 26/07/17)
