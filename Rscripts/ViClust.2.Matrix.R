######################################################################################
###### ViClust results files format to Matrix SNP format ######
###### R script by Pauline Garnier-Géré, June 2015
###### loads ViClust result file and transform it into a classical matrix genotype file format

###### input file= ViClust result file (the most complete one ie res.with.hclust.txt in example):  header == sampleId.new	assayId	ori.genotypeId	Angle	magnitude	call.filter1	HClust.groups
##### output file: header sampleId.new so1.a1 so1.a2 f.so1 shc1.a1 shc1.a2 so2.a1 so2.a2 f.so2 shc2.a1 shc2.a2 
##### with   sampleId.new: ID of individual/well on plate
#####        so1.a1 : allele 1 at original sequenom call for snp 1  (either "A", "G", "C", "T" or "NA" for missing data)
#####        so1.a2 : allele 2 at original sequenom call for snp 1
#####        f.so1: either "call" or "nocall" after filter 1 on original snp 1 call by the sequenom
#####        shc1.a1 : allele 1 based on Ward Hierarchical Clustering call for snp 1
#####        shc1.a2 : allele 2 based on Ward Hierarchical Clustering call for snp 1
#####        so2.a1 : allele 1 at original sequenom call for snp 2  (either "A", "G", "C", "T" , "WC" or "NA" for missing data)
#####        so2.a2 : allele 2 at original sequenom call for snp 2
#####        etc...
#####        Warning: if shcx.ly alleles are proposed as "WC" , it usually means that despite being "NA" based on the sequenom call, they can be grouped with
#####        one particular allelic class based on the Ward clustering method. This is done automatically yet but needs to be modified with the corresponding
#####        allele by the user or be excluded & considered as missing data.

library(seqinr)

#getwd() # check wking dir
## Put the program in the working directory
## Import output file from the Sequenom (ex. with example file provided res.with.hclust.txt)
snpdat<-read.table("res.with.hclust.txt",header=T,dec=".",na.strings="NA",sep="\t",strip.white=T)
snpdat$ori.genotypeId<-as.vector(snpdat$ori.genotypeId)
snpdat$HClust.groups<-as.vector(snpdat$HClust.groups)
for (i in 1:dim(snpdat)[1]) {
    if (is.na(snpdat$ori.genotypeId[i])==FALSE){
    if (snpdat$ori.genotypeId[i]=="A") { snpdat$ori.genotypeId[i]<-"AA"}
    if (snpdat$ori.genotypeId[i]=="G") { snpdat$ori.genotypeId[i]<-"GG"}
    if (snpdat$ori.genotypeId[i]=="C") { snpdat$ori.genotypeId[i]<-"CC"}
    if (snpdat$ori.genotypeId[i]=="T") { snpdat$ori.genotypeId[i]<-"TT"}
    }
    if (is.na(snpdat$HClust.groups[i])==FALSE) {
    if (snpdat$HClust.groups[i]=="A") { snpdat$HClust.groups[i]<-"AA"}
    if (snpdat$HClust.groups[i]=="G") { snpdat$HClust.groups[i]<-"GG"}
    if (snpdat$HClust.groups[i]=="C") { snpdat$HClust.groups[i]<-"CC"}
    if (snpdat$HClust.groups[i]=="T") { snpdat$HClust.groups[i]<-"TT"}
    }
   }

## printing snp list and individuals list
snp.list<-data.frame(seq(1:length(levels(snpdat$assayId))),levels(snpdat$assayId))
colnames(snp.list)<-c("snp.nb","snp.name")
write.table(snp.list,"snp.list.txt", quote=FALSE, row.names=FALSE, sep="\t")
ind.list<-levels(snpdat$sampleId.new)
write.table(as.data.frame(ind.list),"ind.list.txt", quote=FALSE, row.names=FALSE, sep="\t")

## defining structure data for snp as matrix
ind.nb<-length(levels(snpdat$sampleId.new))
snp.nb<-length(levels(snpdat$assayId))

snp.ori.list<-NULL   # initializing colnames for original snp & alleles base call 
snp.loc<-NULL
snp.2all<-NULL             
for (i in 1:snp.nb) {
    snp.loc<-c(snp.loc,paste("so",i,sep="")) 
    for (j in 1:2) {  # always 2 alleles for diploids , whatever the type
    snp.2all<-c(snp.2all,paste(snp.loc,".a",j,sep=""))
    }
   snp.ori.list<-c(snp.ori.list,snp.2all)
   snp.loc<-NULL                          
   snp.2all<-NULL 
  }
  
snp.filter.list<-NULL  # initializing colnames for filter "call"/"nocall" at each snp
for (i in 1:snp.nb) {
    snp.filter.list<-c(snp.filter.list,paste("f.so",i,sep="")) 
    }

snp.hc.list<-NULL   # initializing colnames for snp & alleles Hclust call 
snp.loc<-NULL
snp.2all<-NULL             
for (i in 1:snp.nb) {
    snp.loc<-c(snp.loc,paste("shc",i,sep="")) 
    for (j in 1:2) {
    snp.2all<-c(snp.2all,paste(snp.loc,".a",j,sep=""))
    }
   snp.hc.list<-c(snp.hc.list,snp.2all)
   snp.loc<-NULL                          
   snp.2all<-NULL 
  }

all.col.list<-NULL      # Initializing full matrix column list defined as in initial comments (one snp variables after another)
for (i in 1:snp.nb) {
    all.col.list<-c(all.col.list,snp.ori.list[2*i-1],snp.ori.list[2*i],snp.filter.list[i],snp.hc.list[2*i-1],snp.hc.list[2*i])
    }

## Creating a matrix format structure for each snp with the maximum nb of individuals
snp.mat.all<-matrix(NA,nrow=ind.nb,ncol=0, dimnames=list(ind.list,NULL))
snp.mat.all<-cbind(snp.mat.all,rownames(snp.mat.all))
colnames(snp.mat.all)<-"sampleId.new" 

for (snp in snp.list$snp.name) {  # loop across snps
    snp.i<-snpdat[which(snpdat$assayId==snp),]
    ind.nb.s<-dim(snp.i)[1]
    snp.num<-snp.list[which(snp.list$snp.name==snp),1]
    col.list<-c(paste("s",snp.num,".col1",sep=""),paste("s",snp.num,".col2",sep=""),paste("s",snp.num,".col3",sep=""),paste("s",snp.num,".col4",sep=""),paste("s",snp.num,".col5",sep=""))
    snp.mat<-matrix(NA,nrow=ind.nb.s,ncol=5,dimnames=list(levels(as.factor(as.vector(snp.i$sampleId.new))),col.list)) # create matrix structure with individuals ID 
    snp.mat<-cbind(snp.mat,rownames(snp.mat))
    colnames(snp.mat)[6]<-"sampleId.new" 

    for (mlin in 1:ind.nb.s) { 
      snp.mat[mlin,1] <-unlist(strsplit(as.vector(snp.i$ori.genotypeId[mlin]), NULL))[1]
      snp.mat[mlin,2]<-unlist(strsplit(as.vector(snp.i$ori.genotypeId[mlin]), NULL))[2]
      snp.mat[mlin,3]<-as.vector(snp.i$call.filter1[mlin]) 
      if (is.na(as.vector(snp.i$HClust.groups[mlin]))==FALSE & as.vector(snp.i$HClust.groups[mlin])!="WC") {
      snp.mat[mlin,4]<-unlist(strsplit(as.vector(snp.i$HClust.groups[mlin]), NULL))[1] 
      snp.mat[mlin,5]<-unlist(strsplit(as.vector(snp.i$HClust.groups[mlin]), NULL))[2] 
      }
      if (is.na(as.vector(snp.i$HClust.groups[mlin]))==FALSE & as.vector(snp.i$HClust.groups[mlin])=="WC") {
      snp.mat[mlin,4]<-as.vector(snp.i$HClust.groups[mlin]) 
      snp.mat[mlin,5]<-as.vector(snp.i$HClust.groups[mlin]) 
      }                       
      }
  snp.mat.all<-merge(snp.mat.all,snp.mat,by.x="sampleId.new",by.y="sampleId.new",all.x=T, all.y=T, sort=T)    
    }
head(snp.mat.all)  
colnames(snp.mat.all)[2:dim(snp.mat.all)[2]]<-all.col.list

## extract complete matrix format 
write.table(snp.mat.all,file="snp.mat.all.txt",append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)

## extract matrix format for sequenom base call
write.table(snp.mat.all[,c("sampleId.new",snp.ori.list)],file="snp.mat.ori.txt",append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)


## extract matrix format for Hclust method
write.table(snp.mat.all[,c("sampleId.new",snp.hc.list)],file="snp.mat.hc.txt",append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)


