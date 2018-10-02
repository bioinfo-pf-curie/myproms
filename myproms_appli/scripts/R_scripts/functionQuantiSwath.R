################################################################################
# functionQuantiSwath.R         2.0.3                                          #
# Authors: Matthieu Lhotellier & Alexandre Sta (Institut Curie)                #
# Contact: myproms@curie.fr                                                    #
# Function of quantiSwath.R                                                    #
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



checkdata=function(Data){
  if(sum(colnames(Data)%in%c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity"))!=10){
    stop("Format not good")
  }else{
    return(Data[,c("ProteinName","PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType","Condition","BioReplicate","Run","Intensity")])
  }
}

#################################################################################
## Counting Missing Values
##
## Arguments:
## "areaions" data 
## "step" step of processing (ex: "Before_Normalization")
##
#################################################################################
MissingSample <- function(areaions,step){
  # Build the text file : percentageNAsample.txt
  areaions %>% group_by(Condition,BioReplicate,Run) %>% dplyr::summarise(percentNA = sum(is.na(Intensity))*100/n()) %>%
    unite(Run,c("Condition","BioReplicate","Run"),sep="_") %>% t %>% 
    write.table(.,file.path("results",paste(step,"percentageNAsample.txt",sep="_")),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  
  # Build the image colmissing.jpeg
  areaions %>% group_by(Condition,BioReplicate,Run) %>% dplyr::summarise(Percentage = sum(is.na(Intensity))*100/n()) %>%
    unite(Sample,c("Condition","BioReplicate","Run"),sep="_",remove=FALSE) %>% 
    ggplot(aes( x = Sample , y = Percentage , fill = Condition )) + geom_bar(stat="identity")+ theme(axis.text.x=element_text(angle = 45,vjust=0.5)) + 
    ylim(c(0, 100))
  ggsave(file=file.path("results","graph",paste(step,"colmissing.jpeg",sep="")),width=11,height=9)
  graphics.off()
  
  # Build the image rowmissing.jpeg
  areaions %>% group_by(PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge) %>% dplyr::summarise(nNa=sum(is.na(Intensity))) %>%
    group_by(nNa) %>% dplyr::summarise(n=n()) %>% mutate(Percentage=n*100/sum(n)) %>% 
    ggplot(aes(nNa,Percentage)) + geom_histogram(stat="identity",fill="purple") + ylim(c(0, 100)) + xlab("Number of NA") + theme(legend.position = "none") 
  ggsave(file=file.path("results","graph",paste(step,"nbrowmissing.jpeg",sep="")),width=11,height=9)
  graphics.off()
}


#################################################################################
## Descriptive statistics
##
## Arguments:
## "object" data 
## "step" step of processing (ex: "Before_Normalization")
##
#################################################################################
DataDescPlot=function(object,step){
  BoxplotPlot(object,step)
  DensityPlot(object,step)
  #MAPlot(QuantData,step)
}


#################################################################################
## Boxplot of Sample
##
## Arguments:
## "object" data 
## "step" step of processing (ex: "Before_Normalization")
##
#################################################################################

BoxplotPlot=function(data,step){
  data %>% mutate(Intensity=log2(Intensity)) %>% 
    ggplot(aes(Run,Intensity,fill=Condition)) +
    geom_boxplot(outlier.shape=1,outlier.size=1.5)+
    ggplot2::ylab("Log2-intensities")+xlab("MS runs")
  ggsave(file=file.path("results","graph",paste(step,"allpeptide.jpeg",sep="")),width=11,height=9)
  graphics.off()
}

#################################################################################
## Density of Sample
##
## Arguments:
## "object" data 
## "step" step of processing (ex: "Before_Normalization")
##
#################################################################################


DensityPlot=function(tab,step){
  tab %>% mutate(Intensity=log2(Intensity)) %>% 
    ggplot(aes(Intensity,factor=Run,color=Condition)) + geom_density() + xlab("Log2-Intensity") + ylab("Density")
  ggsave(file=file.path("results","graph",paste(step,"alldensity.jpeg",sep="")),width=11,height=9)
  graphics.off()
}

#################################################################################
## Matrix of correlation of Sample
##
## Arguments:
## "object" data 
## "step" step of processing (ex: "Before_Normalization")
## Supplmentary information :  Agregate on the mean the technical replicates
## 
#################################################################################

MatCorPlot=function(tab,step){
  setDT(tab)
  res=dcast.data.table(tab,PROTEIN+FEATURE~GROUP_ORIGINAL+SUBJECT_ORIGINAL+RUN+LABEL,value.var="ABUNDANCE",fun.aggregate=mean)
  res=setDF(res)
  res=res[,-c(1,2)]
  table=cor(res,use="pairwise.complete.obs",method="pearson")
  write.table(table,file.path("results",paste(step,"matrixcorrelation.txt",sep="")),quote=FALSE,row.names=TRUE,col.names=NA,sep="\t")
  M3=melt(table)
  colnames(M3)<-c("Var1","Var2","R")
  p <- ggplot(M3, aes(Var1, Var2)) + geom_tile(aes(fill =R),colour = "white",limits=c(0,1))
  p <- p+ scale_fill_gradient(low = "white",high = "steelblue",limits=c(0,1))+xlab("") +ylab("")+ scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) + geom_text(aes(fill = R, label = round(R, 2)))
  p + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 9 * 0.8, angle = 270, hjust = 0, colour = "grey50"))
  ggsave(file=file.path("results","graph",paste(step,"matrixcorrelation.jpeg",sep="")),width=11,height=9)
  graphics.off()
  setDF(tab)
}
#################################################################################
## MAplot between condition
##
## Arguments:
## "tab" data 
## "step" step of processing (ex: "Before_Normalization")
##
#################################################################################


MAPlot=function(tab,step){
  setDT(tab)
  res=dcast.data.table(tab,PROTEIN+FEATURE~GROUP_ORIGINAL+SUBJECT_ORIGINAL+RUN,value.var="ABUNDANCE")
  res=setDF(res)
  res=res[,-c(1,2)]
  mva.pairs(na.omit(res),log=FALSE)
   
  
  if(ncol(res)<=9){
    png(filename = paste(path, file.path("results","graph","MAplotPep_A"),step,".png",sep=""),width=1024, height=768)
    mva.pairs(na.omit(res),log=FALSE)
    dev.off()
  }
  setDF(tab)
}

#################################################################################
## PCA transition peptide intensities
##
## Arguments:
## "tab" data 
## "step" step of workflow (ex: "Before_Normalization")
##
#################################################################################

PCAsample=function(tab,step){
  setDT(tab)
  res=dcast.data.table(tab,PROTEIN+FEATURE~GROUP_ORIGINAL+SUBJECT_ORIGINAL+RUN,value.var="ABUNDANCE")
  res=setDF(res)
  res=res[,-c(1,2)]
  pca=PCA(t(na.omit(res)),graph=FALSE,scale.unit = FALSE)
  coord=pca$ind$coord[,1:2]
  colnames(coord)=paste("PC",1:2,sep="")
  scores <- data.frame(rownames(coord), coord[,1:2])
  colnames(scores)[1]="Sample"
  scores[,1]=do.call("rbind",str_split(scores[,1],"_"))[,1]
  p <- ggplot(data=scores, aes(x=PC1, y=PC2,colour=Sample)) + geom_point(size=2.5, alpha=.6) +xlab(paste("PC1 (",round(pca$eig[1,2],2),"%)",sep=""))+ylab(paste("PC2 (",round(pca$eig[2,2],2),"%)",sep=""))+ labs(colour = "Sample")
  p
  ggsave(file=file.path("results","graph",paste(step,"pcasample.jpeg",sep="")),width=11,height=9)
  graphics.off()
  setDF(tab)
}

#################################################################################
## clustering correlation transition peptide intensities
##
## Arguments:
## "tab" data 
## "step" step of workflow (ex: "Before_Normalization")
##
#################################################################################

Clusteringsample=function(tab,step){
  setDT(tab)
  res=dcast.data.table(tab,PROTEIN+FEATURE~GROUP_ORIGINAL+SUBJECT_ORIGINAL+RUN,value.var="ABUNDANCE",function(b)median(b,na.rm=TRUE))
  res=setDF(res)
  res=res[,-c(1,2)]
  #res=na.omit(res)
  dissim=cor(res)
  dissim=1-dissim
  hc=hclust(as.dist(dissim),method="ward.D") # AS
  dd.row<-as.dendrogram(hc)
  ddata_x <- dendro_data(dd.row)
  p2 <- ggplot(segment(ddata_x)) +geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  labs <- label(ddata_x)
  labs$group=as.vector(do.call("rbind",str_split(labs$label,"_"))[,1])
  Sample=labs$group
  ddata_x$labels$label=apply(do.call("rbind",str_split(labs$label,"_"))[,2:3],1,paste,collapse="_")
  ddata_x$labels=data.frame(ddata_x$labels,Sample=labs$group)
  p2<-p2+geom_text(data=label(ddata_x),aes(label=label, x=x, y=0, color=Sample))+xlab("")+ylab("")+ scale_x_discrete(breaks=NULL)
  p2
  ggsave(file=file.path("results","graph",paste(step,"clustersample.jpeg",sep="")),width=11,height=9) 
  graphics.off()
  setDF(tab)
}

#################################################################################
## Analysis Differential and quantification ( modele lineair MSstats)
##
## Arguments:
## "tab" data 
##
#################################################################################


AnalysisDiff=function(tab){
  QuantData=tab
  design<-model.matrix(~0+levels(QuantData$GROUP_ORIGINAL))
  colnames(design)<-unique(levels(QuantData$GROUP_ORIGINAL))
  combF<-combn(colnames(design),2)
  labelcombin<-paste(combF[2,],combF[1,],sep="-")
  Nbcombexpin<-length(labelcombin)
  comparison=c()
  for(i in 1:Nbcombexpin){
    contrast<-paste(combF[2,i],combF[1,i],sep="-")
    contrast.matrix <- as.vector(makeContrasts(contrasts=contrast, levels=design))
    comparison=rbind(comparison,contrast.matrix)
    rownames(comparison)[i]<-labelcombin[i]
  }
  testResultComparison<-groupComparison(contrast.matrix=comparison, data=QuantData, labeled=FALSE, scopeOfBioReplication="restricted", scopeOfTechReplication="expanded", interference=TRUE,featureVar=FALSE,missing.action="nointeraction")  
  return(testResultComparison)
}

AlltopProteins=function(tab,nb){
  return(tab$ComparisonResult[order(tab$ComparisonResult$adj.pvalue),][1:nb,])
}

ComptopProteins=function(tab,nb){
  res=dlply(tab$ComparisonResult,"Label",function(x){
    x[order(x$adj.pvalue),][1:nb,]
  })
  return(do.call("rbind",res))
}

#################################################################################
## Export Data
##
## Arguments:
## "tab" data 
## "step" of workflow pvalue
##
#################################################################################


FormatDataToExport=function(tab,step){
  if(step=="After_correction"){
    res=dlply(tab,"Label",function(x){
      res1=dcast(x,Protein~Label,value.var="log2FC")
      colnames(res1)[2:ncol(res1)]=paste("log2FC_",colnames(res1)[2:ncol(res1)],sep="")
      res2=dcast(x,Protein~Label,value.var="adj.pvalue")
      colnames(res2)[2:ncol(res2)]=paste("pval_",colnames(res2)[2:ncol(res2)],sep="")
      res3=dcast(x,Protein~Label,value.var="SE")
      colnames(res3)[2:ncol(res3)]=paste("SD_",colnames(res3)[2:ncol(res3)],sep="")
      res=cbind(res1,res2[,ncol(res2)],res3[,ncol(res3)])
      colnames(res)[3:ncol(res)]<-c(colnames(res2)[ncol(res2)],colnames(res3)[ncol(res3)])
      return(res)
      })
  }else{
    res=dlply(tab,"Label",function(x){
      res1=dcast(x,Protein~Label,value.var="log2FC")
      colnames(res1)[2:ncol(res1)]=paste("log2FC_",colnames(res1)[2:ncol(res1)],sep="")
      res2=dcast(x,Protein~Label,value.var="pvalue")
      colnames(res2)[2:ncol(res2)]=paste("pval_",colnames(res2)[2:ncol(res2)],sep="")
      res3=dcast(x,Protein~Label,value.var="SE")
      colnames(res3)[2:ncol(res3)]=paste("SD_",colnames(res3)[2:ncol(res3)],sep="")
      res=cbind(res1,res2[,ncol(res2)],res3[,ncol(res3)])
      colnames(res)[3:ncol(res)]<-c(colnames(res2)[ncol(res2)],colnames(res3)[ncol(res3)])
      return(res)
    })
  }
  namevar=llply(res,function(x)return(as.vector(names(x))))
  namevar=unlist(namevar)
  nameprot=llply(res,function(x)return((unique((x$Protein)))))
  nameprot=unique(unlist(nameprot))
  res=llply(res,function(x){
    prot=unique(x[,1])
    naprot=nameprot[which(nameprot%in%prot==FALSE)]
    nb=length(naprot)
    nares=as.data.frame(matrix(NA,nrow=nb,ncol=4))
    colnames(nares)=colnames(x)
    nares[,1]=naprot
    condata=rbind(x,nares)
    condata=condata[order(condata[,1]),]
    return(condata)
  })
  res=do.call("cbind",res)
  colnames(res)=namevar
  if(ncol(res)>4){
    removeprot=grep("Protein",colnames(res))[-1]
    res=(res[,-removeprot])
  }
  return(res)
}

#################################################################################
## Histogramme pvalue
##
## Arguments:
## "res" format export data
## "maxcount" 
## "step" step of worflow
#################################################################################


DescrPvalue=function(res,maxcount,step){
  res=res[,c(1,grep("pval",colnames(res)))]
  res=melt(res,id.var="Protein")
  pdf(file = file.path("results","graph",paste(step,"histpvalue.pdf",sep="")))
  d_ply(res,"variable",function(x){
    m <- ggplot(x, aes(x = value,factor=variable,fill=variable))
    m <- m +geom_histogram(binwidth=0.01, linetype=1, alpha=.25, color="black", fill="peru")+ggtitle(paste("Histogram",unique(x$variable)))+ylab("Frequency")+xlab("P value")+ylim(0,maxcount)
    print(m)
  })
  dev.off()
}

#################################################################################
## Remove Cv outliers transition
##
## Arguments:
## "tab" data
#################################################################################

HandlingCVout=function(tab){
  
  design=unique(tab[,c("GROUP_ORIGINAL","SUBJECT_ORIGINAL","RUN")])
  nbrep=ddply(design,"GROUP_ORIGINAL",function(x)length(x$RUN))
  if(any(nbrep[,2]>1)){
    CV=aggregate(ABUNDANCE~PROTEIN+FEATURE+GROUP_ORIGINAL,data=tab,FUN=function(x)return(sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)))
    colnames(CV)[4]="Value"
    CVout=.outCV(CV)
    res=boxplot(CV[,ncol(CV)],plot=FALSE)
    x=tab[,c("PROTEIN","FEATURE","GROUP_ORIGINAL")]
    setDT(x)
    ID=x[,list(ID=paste(PROTEIN,FEATURE,GROUP_ORIGINAL,sep="_"))]
    rm(x)
    x=CV[which(CV[,ncol(CV)]>res$stats[5,1]),][,c("PROTEIN","FEATURE","GROUP_ORIGINAL")]
    setDT(x)
    IDout=x[,list(ID=paste(PROTEIN,FEATURE,GROUP_ORIGINAL,sep="_"))]
    rm(x)
    setDF(ID)
    setDF(IDout)
    removeCV2=tab[-which(ID[,1]%in%IDout[,1]),]
    return(removeCV2)
  }else{
    return(tab)
  }
}

research=function(x){
  Box=boxplot(x$Value,plot=FALSE)
  M=x[which(x$Value%in%Box$out==TRUE),]
  return(M)
}

.outCV=function(x,path=""){
  CVout<-research(x)
  supCV<-boxplot(x$Value,plot=FALSE)$stats[5,1]
  write.table(data.frame(threshold.CV_Sup=supCV),paste(path,file.path("results","thresholdCV.txt"),sep=""),quote=FALSE,row.names=FALSE,sep="\t")
  
  CV=data.frame(x,id=1:nrow(x))
  png(filename = paste(path,file.path("results","graph","pep_cv_hist.png"),sep=""),width = 800)
  p1 <- ggplot(CV,aes(x=id, y= Value))
  p1 <- p1 + geom_boxplot(color="purple")+ylab("CV")+xlab("") +scale_x_discrete(breaks=NULL)
  p2 <- ggplot(CV,aes(x = Value,factor=factor(0)))
  p2 <- p2 + geom_density(color="purple")+xlab("CV")
  p2<- p2 + theme(legend.position="none") #  AS+ geom_vline(xintercept = c(0,supCV),colour="red")
  grid.arrange(p1,p2,ncol=2, widths=c(1/3,2/3))
  dev.off()
  rm(CV)
  
  write.table(CVout,paste(path,file.path("results","CVoutPep.txt"),sep=""),quote=FALSE,row.names=FALSE,sep="\t")
  return(CVout)
}

#################################################################################
## Grubbs test for peptide (need to do on the levels transition)
##
## Arguments:
## "sub" data
## 
#################################################################################

detectpepout=function(sub){
  sub=na.omit(sub)
  setDT(sub)
  indexprot=sub[,list(indexprot=paste(GROUP_ORIGINAL,RUN,PROTEIN,sep="_"))]
  setDF(sub)
  setDF(indexprot)
  indexprot=as.vector(indexprot)
  colnames(sub)[ncol(sub)]="M"
  sub=data.frame(sub,indexprot=indexprot)
  
  nb=aggregate(M~indexprot,data=sub,FUN=length)
  removeprot=nb[which(nb[,ncol(nb)]<=2),1]
  keepprot=nb[which(nb[,ncol(nb)]>2),1]
  removeprot=sub[which(sub$indexprot%in%removeprot),]
  keepprot=sub[which(sub$indexprot%in%keepprot),]
  
  setDT(keepprot)
  res=keepprot[,list(FEATURE=FEATURE,M=testOut.Grubbs(M,0.05)),by=c("indexprot")]
  indexpep=res[,list(indexpep=paste(indexprot,FEATURE,sep="_"))]
  res=data.frame(res,indexpep=indexpep)
  indexpep=keepprot[,list(indexpep=paste(indexprot,FEATURE,sep="_"))]
  keepprot=data.frame(keepprot,indexpep=indexpep)
  keepprot[which(keepprot[,"indexpep"]%in%res[which(is.na(res$M)),"indexpep"]),"M"]<-NA
  keepprot=keepprot[,-ncol(keepprot),]
  prot=rbind(removeprot,keepprot)
  prot=prot[,-ncol(prot)]
  colnames(prot)[ncol(prot)]="ABUNDANCE"
  return(na.omit(prot))
}

testOut.Grubbs <- function(vect,threshold.out)
{
  if( (length(vect)>2)&((sum(vect[1]==vect))!=sum(!is.na(vect)))){
    p1=grubbs.test(vect,opposite=FALSE,two.sided=FALSE,type=10)$p.value
    p2=grubbs.test(vect,opposite=TRUE,two.sided=FALSE,type=10)$p.value
    h1=grubbs.test(vect,opposite=FALSE,two.sided=FALSE,type=10)$alternative
    alter=paste("highest value", max(vect,na.rm=TRUE), "is an outlier",sep=" ")
    while(((p1<threshold.out) | (p2<threshold.out) )& (length(vect[!is.na(vect)])>3)&sum(vect[!is.na(vect)][1]==vect,na.rm=TRUE)!=length(vect[!is.na(vect)])){
      if(identical(h1,alter)){
        vect[vect==max(vect,na.rm=TRUE)]=NA
        p1=grubbs.test(vect,opposite=FALSE,two.sided=FALSE,type=10)$p.value
        p2=grubbs.test(vect,opposite=TRUE,two.sided=FALSE,type=10)$p.value
        h1=grubbs.test(vect,two.sided=FALSE,type=10)$alternative
        alter=paste("highest value", max(vect,na.rm=TRUE), "is an outlier",sep=" ")
      }
      else{
        vect[vect==min(vect,na.rm=TRUE)]=NA
        p1=grubbs.test(vect,opposite=FALSE,two.sided=FALSE,type=10)$p.value
        p2=grubbs.test(vect,opposite=TRUE,two.sided=FALSE,type=10)$p.value
        h1=grubbs.test(vect,two.sided=FALSE,type=10)$alternative
        alter=paste("highest value", max(vect,na.rm=TRUE), "is an outlier",sep=" ")
      }
    }
  }
  return(vect)
}

#################################################################################
## heatmap pvalue
##
## Arguments:
## "tab" data
## 
#################################################################################


Heatmapprot=function(tab){
  y=dcast(tab,Protein~Label,value.var="adj.pvalue")
  if(ncol(y)==2){
    return(print("no heatmap for one comparison"))
  }
  rownames(y)=y[,1]
  y=y[,-1]
  #rownames(y)=gsub(pattern="sp|","",rownames(y),fixed=TRUE)
  #rownames(y)=gsub(pattern="|","_",rownames(y),fixed=TRUE)
  x=str_split(rownames(y),"_")
  x=llply(x,function(y)return(y[1]))
  x=do.call("rbind",x)
  rownames(y)=x[,1]
  breaks = seq(0,max(y,na.rm=TRUE),length.out=1000)
  gradient1 = colorpanel( sum( breaks[-1]<=0.05 ), "red","orange","orange")
  gradient2 = colorpanel( sum( breaks[-1]>0.05 ), "yellow","yellow","yellow")
  hm.colors = c(gradient1,gradient2)
  pdf(file = file.path("results","graph",paste("heatmap_pvalue.pdf",sep="")),width=6,height=200)
  heatmap.2(as.matrix((y)),dendrogram="none",notecol="black",breaks=breaks,symbreaks=FALSE,col=hm.colors,scale="none",key=TRUE, keysize=1.5,density.info="none", trace="none", cexRow=0.7,cexCol=1.2,main="DE proteines",xlab="Comparison",ylab="Protein",srtCol=45,na.color="black",Colv=FALSE,Rowv=FALSE)
  dev.off()
  
  setDF(tab)
  tab=tab[which(abs(tab$log2FC)>1),]
  y=dcast(tab,Protein~Label,value.var="adj.pvalue")
  rownames(y)=y[,1]
  y=y[,-1]
  rownames(y)=gsub(pattern="sp|","",rownames(y),fixed=TRUE)
  rownames(y)=gsub(pattern="|","_",rownames(y),fixed=TRUE)
  x=str_split(rownames(y),"_")
  x=llply(x,function(y)return(y[1]))
  x=do.call("rbind",x)
  rownames(y)=x[,1]
  breaks = seq(0,max(y,na.rm=TRUE),length.out=1000)
  gradient1 = colorpanel( sum( breaks[-1]<=0.05 ), "red","orange","orange")
  gradient2 = colorpanel( sum( breaks[-1]>0.05 ), "yellow","yellow","yellow")
  hm.colors = c(gradient1,gradient2)
  pdf(file = file.path("results","graph",paste("heatmap_pvalue_cutoff.pdf",sep="")))
  heatmap.2(as.matrix((y)),dendrogram="none",notecol="black",breaks=breaks,symbreaks=FALSE,col=hm.colors,scale="none",key=TRUE, keysize=1.5,density.info="none", trace="none", cexRow=0.7,cexCol=1.2,main="DE proteines",xlab="Comparison",ylab="Protein",srtCol=45,na.color="black",Colv=FALSE,Rowv=FALSE)
  dev.off()
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
## remove protein present in one condition
##
## Arguments:
## "tab" data
## "level" levels protein
## 
#################################################################################

# removeonecond=function(QuantData,level="prot"){
#   QuantData=na.omit(QuantData)
#   setDT(QuantData)
#   if(level=="prot"){
#     nb=QuantData[,list(nb=length(unique(GROUP_ORIGINAL))),by=c("PROTEIN")]
#     setDF(nb)
#     setDF(QuantData)
#     protonecond=nb[which(nb$nb==1),1]
#     if(length(protonecond)>0){
#       index=which(QuantData$PROTEIN%in%protonecond)
#       sub=QuantData[index,]
#       #pdf("Protein_in_one_condition_boxplot.pdf")
#       #d_ply(sub,"PROTEIN",function(y){
#       #  BoxplotPlotSub(y)
#       #})
#       #dev.off()
#       BoxplotPlot(sub,"Protein_in_one_condition_After_out")
#       write.table(sub,file.path("results",paste("ProteinsOneCondition.txt",sep="")),quote=FALSE,row.names=FALSE,sep="\t")
#       QuantData=QuantData[-which(QuantData$PROTEIN%in%protonecond),]
#     }
#   }else{
#     tab[,list(nb=length(unique(PEPTIDE))),by=c("PROTEIN","GROUP_ORIGINAL")]
#     nb=QuantData[,list(nb=length(unique(PEPTIDE))),by=c("PROTEIN","GROUP_ORIGINAL")]
#   }
#   return(QuantData)
# }

#################################################################################
## Boxplot
##
## Arguments:
## "tab" data
## 
#################################################################################

BoxplotPlotSub=function(QuantData){
  QuantData=QuantData[which(QuantData$PROTEIN==unique(QuantData$PROTEIN)),]
  QuantData$ABUNDANCE=as.numeric(QuantData$ABUNDANCE)
  QuantData$RUN=as.factor(as.numeric(as.vector(QuantData$RUN)))
  m <- ggplot(QuantData,aes(y = ABUNDANCE,x=RUN))
  tempGroupName<-unique(QuantData[,c("GROUP_ORIGINAL","RUN")])
  groupAxis<-as.numeric(xtabs(~GROUP_ORIGINAL,tempGroupName))
  cumGroupAxis<-cumsum(groupAxis)
  lineNameAxis<-cumGroupAxis[-nlevels(QuantData$GROUP_ORIGINAL)]
  groupName<-data.frame(RUN=c(0,lineNameAxis)+groupAxis/2+0.5,y=rep(30-1,length(groupAxis)),Name=levels(QuantData$GROUP_ORIGINAL))
  m<-m+geom_boxplot(aes_string(fill='SUBJECT_NESTED'),outlier.shape=1,outlier.size=1.5)+annotate("text",x=groupName$RUN,y=groupName$y,label=groupName$Name,size=4,angle=0)+ylab("Log2-intensities")+xlab("MS runs")+ggtitle(paste("Boxplot Protein: ",unique(QuantData$PROTEIN) ,sep=""))+geom_vline(xintercept=lineNameAxis+0.5,colour="grey",linetype="longdash")+ guides(fill=FALSE)+theme_bw()
  plot(m)
}


#################################################################################
## Deviation sample
##
## Arguments:
## "tab" data
## "step" of worklow
## 
#################################################################################

deviationsample=function(tab,step){
  setDT(tab)
  sdintra=tab[,list(intra=sd(ABUNDANCE,na.rm=TRUE)),by=c("PROTEIN","FEATURE","GROUP_ORIGINAL")]
  sdinter=tab[,list(inter=sd(ABUNDANCE,na.rm=TRUE)),by=c("PROTEIN","FEATURE")]
  
  m<-ggplot(sdintra,aes(y=intra,x=GROUP_ORIGINAL,fill=GROUP_ORIGINAL))
  m+geom_boxplot(outlier.shape=1,outlier.size=1.5)+ ylab("SD")+xlab("CONDITION")+guides(fill=FALSE)
  ggsave(file=file.path("results","graph",paste(step,"_deviation_intrasample.jpeg",sep="")),width=11,height=9)
  graphics.off()
  
  qplot(x="ALL",y=inter,data=sdinter,geom="boxplot",fill="ALL")+xlab("")+ylab("SD")+guides(fill=FALSE)
  ggsave(file=file.path("results","graph",paste(step,"_deviation_intersample.jpeg",sep="")),width=11,height=9)
  graphics.off()  
}
#################################################################################
## PCA transition peptide
##
## Arguments:
## "tab" data
## "step" of worklow
## 
#################################################################################

PCAind=function(tab,step){
  setDT(tab)
  res=dcast.data.table(tab,PROTEIN+FEATURE~GROUP_ORIGINAL+SUBJECT_ORIGINAL+RUN,value.var="ABUNDANCE")
  ID=res[,list(ID=paste(PROTEIN,FEATURE,sep="_"))]
  res=setDF(res)
  ID=setDF(ID)
  rownames(res)=ID[,1]
  res=res[,-c(1,2)]
  pca=PCA(na.omit(res),graph=FALSE,scale.unit = FALSE)
  coord=pca$ind$coord[,1:2]
  colnames(coord)=paste("PC",1:2,sep="")
  scores <- data.frame(rownames(coord), coord[,1:2])
  colnames(scores)[1]="Source"
  
  scores=data.frame(scores,orig="NOHUMAN")
  scores$orig="NOHUMAN"
  scores$orig=as.factor(scores$orig)
  scores$orig=as.vector(scores$orig)
  scores[grep("UPS",scores[,1]),"orig"]<-"HUMAN"
  scores$orig=as.factor(scores$orig)
  p <- ggplot(data=scores, aes(x=PC1, y=PC2,colour=orig)) + geom_point(size=2.5, alpha=.6) +xlab(paste("PC1 (",round(pca$eig[1,2],2),"%)",sep=""))+ylab(paste("PC2 (",round(pca$eig[2,2],2),"%)",sep=""))+ labs(colour = "Sample")
  p
  ggsave(file=file.path("results","graph",paste(step,"pcaind.jpeg",sep="")),width=11,height=9)
  graphics.off()
  setDF(tab)
}

#################################################################################
## heatmap log2FC
##
## Arguments:
## "tab" data
## 
#################################################################################


Heatmapratioprot=function(tab){
  setDF(tab)
  tab=dcast(tab,Protein~Label,value.var="log2FC")
  if(ncol(tab)==2){
    return(print("no heatmap for one comparison"))
  }
  tab=melt(tab,id.var=1)
  colnames(tab)[2:3]=c("Label","log2FC")
  tab[which(is.na(tab$log2FC)),"log2FC"]=0
  y=dcast(tab,Protein~Label,value.var="log2FC")
  y=y[,-1]
  breakup=max(abs(ceiling(min(y,na.rm=TRUE)-1)),abs(ceiling(max(y,na.rm=TRUE)+1)))
  breakdown=-breakup
  breaks=breakdown:breakup
  p<- ggplot(tab, aes(Label,Protein))
  p<- p + geom_tile(aes(fill=log2FC), colour="white")
  p<- p + scale_fill_gradientn(colours = c("steelblue", "white", "red" ), breaks=breaks, labels=format(breaks))
  blue.bold.italic.16.text <- element_text(face = "bold.italic", color = "black", size = 4)
  p + theme(axis.text = blue.bold.italic.16.text)+ theme(axis.ticks = element_line(size = 0.01))
  ggsave("results/graph/heatmap_log2FC.pdf",width=4,height=49)
  graphics.off()
}

#################################################################################
## acp log2FC
##
## Arguments:
## "tab" data
## 
#################################################################################


ACPratioprot=function(tab){
  setDT(tab)
  res=dcast.data.table(tab,Protein~Label,value.var="log2FC")
  setDF(res)
  if(ncol(res)==2){
    return(print("no heatmap for one comparison"))
  }
  rownames(res)=res[,1]
  res=res[,-1]
  pca=PCA(na.omit(res),graph=FALSE,scale.unit = FALSE)
  coord=pca$ind$coord[,1:2]
  colnames(coord)=paste("PC",1:2,sep="")
  scores <- data.frame(rownames(coord), coord[,1:2])
  colnames(scores)[1]="Source"
  scores=data.frame(scores,orig="NOHUMAN")
  scores$orig="NOHUMAN"
  scores$orig=as.factor(scores$orig)
  scores$orig=as.vector(scores$orig)
  scores[grep("UPS",scores[,1]),"orig"]<-"HUMAN"
  scores$orig=as.factor(scores$orig)
  p <- ggplot(data=scores, aes(x=PC1, y=PC2,colour=orig)) + geom_point(size=2.5, alpha=.6) +xlab(paste("PC1 (",round(pca$eig[1,2],2),"%)",sep=""))+ylab(paste("PC2 (",round(pca$eig[2,2],2),"%)",sep=""))+ labs(colour = "Sample")
  p
  ggsave(file=file.path("results","graph","pcalog2FC.jpeg"),width=11,height=9)
  graphics.off()
  setDF(tab)
}

#################################################################################
## export pep used for quantification
##
## Arguments:
## "QuantData" data
##
## Description : 
##   List all the peptides(*charge) in each protein*condition for each GROUP_ORIGINAL
##   and RUN. ie, a peptide can be seen many times in many biological or technical
##   sample but we count it once per sample (we does not count the TRANSITION). 
##
#################################################################################


exportpep=function(QuantData){
  setDT(QuantData)
  QuantData$PEPTIDE=gsub(pattern="\\*",replacement=":",x=as.character(QuantData$PEPTIDE)) # AS : cange the "*" to ":" like in the initial data 

  setDF(QuantData)
  tmp = unique(QuantData[,c("PROTEIN","PEPTIDE","GROUP_ORIGINAL","SUBJECT_ORIGINAL","RUN")])
  setDT(tmp)
  usedpep=tmp[,list(ID=paste(PEPTIDE,collapse="|")),by=c("PROTEIN","GROUP_ORIGINAL")]
  # AS :  setDF(QuantData)  
  write.table(usedpep,"results/list_of_peptide_for_quantification.txt",quote=FALSE,row.names=FALSE,sep="\t")
}


#################################################################################
## stard ad end each step
##
## Arguments:
## "QuantData" data
## 
#################################################################################
start <- function(message)
{
 cat(paste("Start : ",message, "\n \n",sep=" "))
 step.begin <- c(message,proc.time()[3])
 return(step.begin)
}

end <- function(step.begin)
{
 cat(paste("End : ",step.begin[1],"\n",sep=" "))
t2 <- proc.time() # Print message for .rout
cat(paste("Time of the step : ",round((t2[3]-as.numeric(step.begin[2])),2),"s\n",sep=""))
cat(paste("Time since the begining : ",round((t2[3]-time.start),2),"s\n",sep=""))
cat("############################# \n \n")
 return()
}

#################################################################################
## select peptide present in each condition
##
## Arguments:
## "QuantData" data
## 
#################################################################################

pepallcondition=function(QuantData){
  setDT(QuantData)
  nb=QuantData[,list(nb=length(unique(GROUP_ORIGINAL))),by=c("PROTEIN","PEPTIDE")]
  nbtable=table(nb[,3,with=FALSE])
  nbtable=(nbtable/sum(nbtable))*100
  setDF(nb)
  nb=(nb[which(nb[,3]==6),])
  setDT(nb)
  IDnb=nb[,list(ID=paste(PROTEIN,PEPTIDE,sep="_"))]
  setDF(IDnb)
  ID=QuantData[,list(ID=paste(PROTEIN,PEPTIDE,sep="_"))]
  setDF(ID)
  setDF(QuantData)
  QuantData=QuantData[which(ID[,1]%in%IDnb[,1]),]
  return(QuantData)
}

pepallcondition2=function(areaions){
  setDT(areaions)
  nbcond=length(unique(areaions$Condition))
  nb=areaions[,list(nb=nbcond),by=c("ProteinName","Peptihead(QuantDdeSequence","PrecursorCharge")]
  nbtable=table(nb[,4,with=FALSE])
  nbtable=(nbtable/sum(nbtable))*100
  setDF(nb)
  nb=(nb[which(nb[,4]==nbcond),])
  setDT(nb)
  IDnb=nb[,list(ID=paste(ProteinName,PeptideSequence,PrecursorCharge,sep="_"))]
  setDF(IDnb)
  ID=areaions[,list(ID=paste(ProteinName,PeptideSequence,PrecursorCharge,sep="_"))]
  setDF(ID)
  setDF(areaions)
  areaions=areaions[which(ID[,1]%in%IDnb[,1]),]
  return(areaions)
}

####>Revision history<####
# 2.0.3 modification of the others functions of visualisation (AS 29/01/18) 
# 2.0.2 modfiication of MissingSample to fit the areaions shape (AS 24/01/18) 
# 2.0.1 correction of the name of the labels in a graph (AS 16/01/18) 
# 2.0.0 the code is now compatible with MSstats >3 (AS 21/11/17) 
# 1.0.7 Fix the bug log2FC, constraint the name of this column to be log2FC and not logFC (AS 29/11/16) 
# 1.0.6 add a function to print start and end of each step (AS 04/08/16)
# 1.0.5 remove : prot from the names of peptides (AS 26/07/16)
# 1.0.4 remove 01 (AS 18/07/16)
# 1.0.3 replace "_" by "|" in "list_of_peptide_for_quantification.txt"  (AS 18/07/16)
# 1.0.2 correction of graphical outputs (AS 06/06/16)
# 1.0.1 correction of graphical outputs (AS 04/11/15)





