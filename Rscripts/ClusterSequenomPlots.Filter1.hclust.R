###### How to cite ViClust:
# Bouteiller XP, Barraquand F, Garnier-Géré P, Harmand N, Laizet Y, Raimbault A, Segura R, Lassois L, Monty A, Verdu C, Mariette S & Porté AJ (2018)
# No evidence for genetic differentiation in juvenile traits between Belgian and French populations of the invasive tree Robinia pseudoacacia .
# Plant Ecology and Evolution 151 (1): 5–17
## See also Website https://github.com/garniergere/ViClust/ for detailed description,  examples and Galaxy version

#############################################################################
# Standalone R script
# deals with raw xml files from Sequenom SNP genotypes
# allows batch visualization and alternative clustering
# Contact: Pauline Garnier-Géré, pauline.garnier-gere@inra.fr
# March 2018
#############################################################################
## Code below is the version for using with the Rscript.exe in command line under windows or linux

## Examples of command lines (NB since Heterozygotes are in black, better not to use the black color for homozygotes clusters
## unless code is changed below accordingly)
##NB: below if command line ran under dos prompt in windows then "" are used if there are space in folders within the used path

## example 1 using a selection of SNPs in SNPlist.ex.txt, all individuals, adjusted scale for each plot (first "no"), Filter 1 at 8 and Hclust plots (before last arg is "yes")
# "C:\Path ToR.exe.folder\Rscript.exe" C:\path.to.working.folder\ClusterSequenomPlots.Filter1.hclust.R viclust.ex.xml SNPlist.ex.txt all red green blue orange no sequenom.plots out filter1.plots 8 yes hclust.plots

## example 2 using all SNPs, all individuals, same scale plots (first yes), no Filter 1 and no Hclust plots (last "no"),  other choice of colors
# "C:\Path ToR.exe.folder\Rscript.exe" C:\path.to.working.folder\ClusterSequenomPlots.Filter1.hclust.R viclust.ex.xml all all red cyan green blue yes sequen.plot2 out.files2 filter1.plot2 0 no

##########################################################################
#rm(list=ls())
library(XML)
##########################################################################
cmd_args <- commandArgs(trailingOnly = TRUE)

## Defining arguments if using Rscript.exe
xmlfile<-cmd_args[1] # argument 1 is the xml file to import
snplist<-cmd_args[2] # argument 2 is the chosen SNP list file to plot, put "all" in command line if all SNP are being analyzed/plotted
indlist<-cmd_args[3] # argument 3 is the individual list file to plot, put "all" if all individuals are being analyzed/plotted
colclusterA<-cmd_args[4] # colour cluster base "A"    red
colclusterG<-cmd_args[5] # colour cluster base "G"    green
colclusterC<-cmd_args[6] # colour cluster base "C"    blue
colclusterT<-cmd_args[7] # colour cluster base "T"    orange
samescale<-cmd_args[8]   # "yes" or "no" for same scale all plots or not, no by default
SequenomPlots <-cmd_args[9] # folder for sequenom plots
output.snpcall<-cmd_args[10] # folder used in Rscript for cluster call results for selected SNPs and individuals
filter1Plots<-cmd_args[11]   #  folder for plots after filter 1
filter1val<-cmd_args[12]     # Minimal magnitude threshold value based on SNR , excluding data with too low signals, 
                             # will plot only called data, default is 6, alternative no or more stringent value  
hclust.condition<-cmd_args[13]   # "yes" or "no" for performing Hclust method, "no" by default
HclustPlots<-cmd_args[14]    #  folder for plots after Hclust (Hierarchical Clustering), if "no" above, can be empty in command line
## create folders (Rscript)
dir.create(SequenomPlots)
dir.create(output.snpcall)
dir.create(filter1Plots)
dir.create(HclustPlots)

###########################################################################
## Loading xml output/result file
###########################################################################
data.xml <- xmlParse(xmlfile)   # arg 1 in Rscript command line
data.list<-xmlToList(data.xml)
records<-t(as.data.frame(data.list$typeranalyzer$'all-records')) ## OK after changing "-" into "." in all-records worked  or put all-records in ''
record2<-records[,c("sampleId","assayId","genotypeId","height","area","mass","peakScore","peakType","snr","resolution","description","wellPosition")]

##########################################################################
### defining list of SNPs and samples (=individuals)
##########################################################################

if (snplist=="all"){  # snplist is argument 2
 SNPlist<-levels(as.factor(record2[,"assayId"]))
}
if (snplist!="all"){ # if snplist argument different from "all", then SNPlist is filled with the list of SNP in the text file given in the snplist argument
 SNPlist<-scan(snplist,sep="", strip.white=TRUE, what="") # for use with Rscript
}
if (indlist=="all"){
INDlist<-levels(as.factor(record2[,"sampleId"]))    # will need to define INDlist2 with sampleId.new when treating pb with replicate individuals 
}
if (indlist!="all"){
 INDlist<-scan(indlist,sep="", strip.white=TRUE, what="")
}

#############################################
### Extract data for chosen SNPs & individuals

selsnp<-matrix(NA,nrow=0,ncol=dim(record2)[2])
colnames(selsnp)<-colnames(record2) 
for (snp in SNPlist) {
     if (snplist=="all")  {selsnp<-record2}
     if (snplist!="all") {
           selsnp.s<-record2[which(record2[,"assayId"]==snp),]    #selsnp.s<-record2[which(record2[,"assayId"]==snp & record2[,"ori.genotypeId"]=="G"),]
           selsnp<-rbind(selsnp,selsnp.s)
     }
}
        
selsnpf<-matrix(NA,nrow=0,ncol=dim(selsnp)[2])
colnames(selsnpf)<-colnames(selsnp) 
for (samp in INDlist)  {
     if (indlist=="all") {selsnpf<-selsnp}
     if (indlist!="all") {
           selsnp.i<-selsnp[which(selsnp[,"sampleId"]==samp),]
           selsnpf<-rbind(selsnpf,selsnp.i)
          }
        } 
selsnp<-selsnpf

### replace genotypeId by ori.genotypeId
colnames(selsnp)[3]<-"ori.genotypeId"

################################################################################
### export raw output from parsed xml file & selected SNP and individual lists
################################################################################
# replace all "" ori.genotypeId by "NA" for export to text file without quote
selsnp.export<-selsnp
for (i in 1:dim(selsnp.export)[1]){
 if (selsnp.export[i,"ori.genotypeId"]=="") {selsnp.export[i,"ori.genotypeId"]<-"NA"} 
 }
write.table(as.data.frame(selsnp.export),file=paste(output.snpcall,"/","rawdat.txt",sep=""),append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)

#######################################
### define sampleID.new in selsnp matrix for dealing with individuals which are replicated
### so that only 3 record by data exist when computing angle & magnitude
### sampleID.new is the merging of the sampleID and of the well number

sampleId.new<-matrix(NA,ncol=1,nrow=dim(selsnp)[1])
colnames(sampleId.new)<-"sampleId.new"
for (i in 1:dim(selsnp)[1]) {
sampleId.new[i,1]<-paste(selsnp[i,"sampleId"],".",selsnp[i,"wellPosition"],sep="")
}
selsnp1<-selsnp
selsnp<-cbind(sampleId.new,selsnp1)
#head(selsnp)
#############################################
## create Matrix with 1 line per combination SNP*individual*well (sampleId.new, see above)
## & computing angle value, magnitude value + NA instead of ""
#############################################

snpres<-matrix(NA,nrow=0,ncol=5)
colnames(snpres)<-c(colnames(selsnp)[c(1,3,4)],"Angle","magnitude") 
# defining complete sample list of non redundant individuals
samplelist<-levels(as.factor(selsnp[,"sampleId.new"]))

 for (snp in SNPlist) {
     for (sam in samplelist) {
      samp<-selsnp[which(selsnp[,"assayId"]==snp & selsnp[,"sampleId.new"]==sam),]  # extract 3 lines for one SNP & one individual (sample)
      if (length(samp)!=0) {   # for non-empty snpres files
         highmass.record<-samp[which(samp[,"mass"]==max(samp[,"mass"])),]  # str(highmass.record)  --> still a matrix structure with colnames
         samp.less.probe<-samp[-which(samp[,"mass"]==min(samp[,"mass"])),]
         lowmass.record<-samp.less.probe[which(samp.less.probe[,"mass"]==min(samp.less.probe[,"mass"])),]
         if (as.numeric(lowmass.record["height"])==0) {
             samp.res<-t(matrix(c(sam,snp,samp[1,"ori.genotypeId"],pi/2,
             sqrt((as.numeric(highmass.record["snr"]))^2 + (as.numeric(lowmass.record["snr"]))^2)),dimnames=list(colnames(snpres))))
         }
         if (as.numeric(lowmass.record["height"])!=0) {
             samp.res<-t(matrix(c(sam,snp,samp[1,"ori.genotypeId"],atan(as.numeric(highmass.record["height"])/as.numeric(lowmass.record["height"])), # use height for angle
                  sqrt((as.numeric(highmass.record["snr"]))^2 + (as.numeric(lowmass.record["snr"]))^2)),dimnames=list(colnames(snpres)))) # but snr for magnitude
         }
         if (samp.res[,"ori.genotypeId"]=="") {samp.res[,"ori.genotypeId"]<-"NA"}
             snpres<-rbind(snpres,samp.res)
      }   # end loop for non empty snpres file
     }   # end loop samples with new ID (so no replicate ID)
 } # end loop across SNPlist

############################################################################################################
## Part 1 plotting SNPs (& assigning them to clusters) as in the Sequenom sofware after the Cluster Call
############################################################################################################

for (snp in SNPlist){  # loop across SNPs
tiff(filename=paste(SequenomPlots,"/",snp,".tiff",sep=""),width=690,height=690)  # SequenomPlots is argument 9, activated when using Rscript
                                                                                 # use "png" instead of "tiff" for less memory use
   snpres.1<-as.data.frame(snpres[which(snpres[,"assayId"]==snp),])
  if (length(snpres.1$sampleId.new)<2) {   # if only 1 datapoint, no Clustering
        plot(0,0,col="white",pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("less than 2 called datapoints for"," ",snp,sep=""))
  }
  if (length(snpres.1$sampleId.new)>=2) {  # if at least 2 called datapoints
       genoclass<-levels(as.factor(as.vector(snpres.1[,"ori.genotypeId"])))
       if (length(genoclass)==0) {   # if no data
         plot(0,0,col="white",pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("No data for"," ",snp,sep=""))
       }

    if (length(genoclass)!=0) {   # start loop at least one cluster  (including No call cluster)
        if (length(which(genoclass=="NA"))==0) {genoclass.nona<-genoclass}
        if (length(which(genoclass=="NA"))!=0) {genoclass.nona<-genoclass[-which(genoclass=="NA")]
        } # remove "NA" levels if they do exist

        ### Plots different group of genotypes with colours as argument in command line
        if (samescale=="no") {
           mxlim<-max(as.numeric(as.vector(snpres.1$Angle))+0.1)
           mylim<-max(as.numeric(as.vector(snpres.1$magnitude))+3)
        }
        if (samescale=="yes") {
           snpres.df<-as.data.frame(snpres)
           mxlim<-max(as.numeric(as.vector(snpres.df$Angle))+0.1)
           mylim<-max(as.numeric(as.vector(snpres.df$magnitude))+3)
        }
        geno.gp<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[1]),]  # define group of points for first geno class
        geno.gp$Angle<-as.numeric(as.vector(geno.gp$Angle))
        geno.gp$magnitude<-as.numeric(as.vector(geno.gp$magnitude))

        cluster.list<-c("A","G","C","T","AG","GA","AC","CA","AT","TA","GC","CG","GT","TG","CT","TC")
        #color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"purple","purple","purple2","purple2","purple4","purple4","cyan","cyan","turquoise2","turquoise2","turquoise4","turquoise4")
        color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"black","black","black","black","black","black","black","black","black","black","black","black") # new color list where all heteroz are black   01/04/15

        for (i in 1:length(cluster.list)) {
          if  (dim(geno.gp)[1]==0) {col1<-"white"}      # is.data.frame(snpres.1) #TRUE
          else if (genoclass.nona[1]==cluster.list[i]) {col1<-color.list[i]}
        }
      # colors()
        if (dim(geno.gp)[1]==0) { # if only NA data, geno.gp has a dim of 0
         plot(0,0,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,mylim),xlim=c(0,mxlim),main=paste("Only NA data for"," ",snp,sep=""))
        # defining & plotting eventual nocall cluster
         geno.NA<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
         geno.NA$Angle<-as.numeric(as.vector(geno.NA$Angle))
         geno.NA$magnitude<-as.numeric(as.vector(geno.NA$magnitude))
         points(geno.NA$Angle,geno.NA$magnitude,col="grey",pch=16,cex=0.8)
        }
     
      if (dim(geno.gp)[1]!=0) { #if at least one cluster which is not NA
        plot(geno.gp$Angle,geno.gp$magnitude,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(-0.1*mylim,mylim),xlim=c(-0.15,mxlim),main=snp)
        legend(-0.1,-1 , legend = c(genoclass.nona[1], "NA"), col = c(col1,"grey"), pch=19, cex=0.9, bty="n", xpd=TRUE, horiz=TRUE)
        # plotting eventual NA datapoints
        geno.NA<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
        geno.NA$Angle<-as.numeric(as.vector(geno.NA$Angle))
        geno.NA$magnitude<-as.numeric(as.vector(geno.NA$magnitude))
        points(geno.NA$Angle,geno.NA$magnitude,col="grey",pch=16,cex=0.8)

       if (length(genoclass.nona)>1) { # if > 1 class with called geno
        geno.gp2<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[2]),]  # defining group 2 & 3 at maximum, can be empty  the plot won't work
        geno.gp2$Angle<-as.numeric(as.vector(geno.gp2$Angle))
        geno.gp2$magnitude<-as.numeric(as.vector(geno.gp2$magnitude))

        for (i in 1:length(cluster.list)) {
          if  (dim(geno.gp2)[1]==0) {col2<-"white"}
          if (genoclass.nona[2]==cluster.list[i]) {col2<-color.list[i]}
        }
        points(geno.gp2$Angle,geno.gp2$magnitude,col=col2,pch=16,cex=0.8)
        legend(-0.1,-1, legend=c("                              "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)  # draws a white rectangle to cover previous legend
        legend(-0.1,-1 , legend = c(genoclass.nona[1],genoclass.nona[2],"NA"), col = c(col1, col2, "grey") , pch = 19,cex = 0.9, bty="n", xpd=TRUE, horiz=TRUE)

        if (length(genoclass.nona)>2) { # case with 3 nona clusters
         geno.gp3<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[3]),]
         geno.gp3$Angle<-as.numeric(as.vector(geno.gp3$Angle))
         geno.gp3$magnitude<-as.numeric(as.vector(geno.gp3$magnitude))
          for (i in 1:length(cluster.list)) {
           if  (dim(geno.gp3)[1]==0) {col3<-"white"}
           if (genoclass.nona[3]==cluster.list[i]) {col3<-color.list[i]}
          }
         points(geno.gp3$Angle,geno.gp3$magnitude,col=col3,pch=16,cex=0.8)
         legend(-0.1,-1, legend=c("                             "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)
         legend(-0.1,-1, legend = c(genoclass.nona[1],genoclass.nona[2],genoclass.nona[3],"NA"), col = c(col1, col2, col3,"grey") , pch = 19,cex = 0.9,bty="n", xpd=TRUE, horiz=TRUE)
        } # end if 3 clusters of genotypes with or without NA
       } # end > 1 cluster call of genotypes to plot with or without NA

    }  # end at least 1 called cluster with or without NA points
   }  # end at least 1 cluster (including NA)
  }  # end at least 2 called datapoints

dev.off()     # closing plot device if tiff fun used above & if using Rscript in command line
}          # end loop across SNPs & SEQUENOM clusters

####################################################################################################
#### part 2 for filtering based on magnitude and plotting data after excluding filtered points
####################################################################################################

#### Create a call factor for excluding datapoints with too low magnitude based on SNR (Filter 1)
  call.filter1<-matrix(NA,ncol=1,nrow=dim(snpres)[1])
  colnames(call.filter1)<-"call.filter1"
## Exporting only Angle & Magnitude
 snpres.df<-as.data.frame(snpres)
 snpres.df$magnitude<-as.numeric(as.vector(snpres.df$magnitude))     
filter1val<-as.numeric(filter1val)
## Condition for no additional filter so use default value of 6
if (filter1val==0) {
  for (i in 1:dim(snpres.df)[1]) {
  if (snpres.df$magnitude[i]<=6) {call.filter1[i,1]<-"nocall.below6" }
  else if (snpres.df$magnitude[i]>6) {call.filter1[i,1]<-"call.above6" }
  }   
  snpres.df2<-cbind(snpres.df,call.filter1)       
  snpres2<-snpres.df2
  write.table(snpres2,file=paste(output.snpcall,"/","res.ang.mag.txt",sep=""),append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)
  }
## condition for additional filter
if (filter1val>0) {
 for (i in 1:dim(snpres.df)[1]) {
  if (snpres.df$magnitude[i]<=filter1val) {call.filter1[i,1]<-"nocall" }
  else if (snpres.df$magnitude[i]>filter1val) {call.filter1[i,1]<-"call" }
 }
 snpres.df2<-cbind(snpres.df,call.filter1)       
 snpres2<-snpres.df2
 write.table(snpres2,file=paste(output.snpcall,"/","res.filter1.txt",sep=""),append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)
}

################################################################
#### Plots after applying filter1 based on SNR #################

if (filter1val==0) { #start condition filter 1
    snpres.f1<-subset(snpres2,subset=call.filter1=="call.above6")  # Df with only "call.above6" points = default
    } # end condition filter 1
if (filter1val>0) {
    snpres.f1<-subset(snpres2,subset=call.filter1=="call")  # Df with only "call" points based on first filter on magnitude
}
# snpres.f1 is the data file after filter 1 based on magnitude

for (snp in SNPlist){  # loop across SNPs with only "call" ori.genotypeId & some NA might still there (grey points)
 tiff(filename=paste(filter1Plots,"/",snp,".tiff",sep=""),width=690,height=690)
       snpres.1<-as.data.frame(snpres.f1[which(snpres.f1[,"assayId"]==snp),])
 if (length(snpres.1$sampleId.new)<2) {   # if only 1 sample, no clustering
    plot(0,0,col="white",pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("less than 2 called datapoints for"," ",snp,sep=""))
 }
 if (length(snpres.1$sampleId.new)>=2) {  # if at least 2 called datapoints
    genoclass<-levels(as.factor(as.vector(snpres.1[,"ori.genotypeId"])))
    if (length(genoclass)==0) {   # if no data
        plot(0,0,col="white",pch=16,cex=1, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("No data for"," ",snp,sep=""))
    }
  if (length(genoclass)!=0) {  # if at least one cluster with maybe NA
    if (length(which(genoclass=="NA"))==0) {genoclass.nona<-genoclass}
    if (length(which(genoclass=="NA"))!=0) {genoclass.nona<-genoclass[-which(genoclass=="NA")]}
    ### Plots different group of genotypes with colors as argument in command line
    if (samescale=="no") {
    mxlim<-max(as.numeric(as.vector(snpres.1$Angle))+0.1)
    mylim<-max(as.numeric(as.vector(snpres.1$magnitude))+3) }
    if (samescale=="yes") {
    snpres.df<-as.data.frame(snpres)
    mxlim<-max(as.numeric(as.vector(snpres.df$Angle))+0.1)
    mylim<-max(as.numeric(as.vector(snpres.df$magnitude))+3) }

    geno.gp<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[1]),]
    geno.gp$Angle<-as.numeric(as.vector(geno.gp$Angle))
    geno.gp$magnitude<-as.numeric(as.vector(geno.gp$magnitude))

    # defining colors for clusters
    cluster.list<-c("A","G","C","T","AG","GA","AC","CA","AT","TA","GC","CG","GT","TG","CT","TC")
    color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"black","black","black","black","black","black","black","black","black","black","black","black") # new color list where all heteroz are black   01/04/15

    for (i in 1:length(cluster.list)) {
      if  (dim(geno.gp)[1]==0) {col1<-"white"}
      else if (genoclass.nona[1]==cluster.list[i]) {col1<-color.list[i]}
    }
    if (dim(geno.gp)[1]==0) { # if only NA data
     plot(0,0,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,mylim),xlim=c(0,mxlim),main=paste("Only NoCall data for"," ",snp,sep=""))
     #  defining & plotting eventual nocall cluster
     geno.nocall<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
     geno.nocall$Angle<-as.numeric(as.vector(geno.nocall$Angle))
     geno.nocall$magnitude<-as.numeric(as.vector(geno.nocall$magnitude))
     points(geno.nocall$Angle,geno.nocall$magnitude,col="grey",pch=16,cex=0.8)
    }

   if (dim(geno.gp)[1]!=0) { #if at least one class which is not NA
    plot(geno.gp$Angle,geno.gp$magnitude,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(-0.1*mylim,mylim),xlim=c(-0.15,mxlim),main=snp)
    legend(-0.1,-1 , legend = c(genoclass.nona[1], "NA"), col = c(col1,"grey"), pch=19, cex=0.9, bty="n", xpd=TRUE, horiz=TRUE) # modif 30/03/15
    # defining & plotting eventual NA cluster
    geno.NA<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
    geno.NA$Angle<-as.numeric(as.vector(geno.NA$Angle))
    geno.NA$magnitude<-as.numeric(as.vector(geno.NA$magnitude))
    points(geno.NA$Angle,geno.NA$magnitude,col="grey",pch=16,cex=0.8)

    if (length(genoclass.nona)>1) { # if more than one class with no NA genotypes
     geno.gp2<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[2]),]   
     geno.gp2$Angle<-as.numeric(as.vector(geno.gp2$Angle))   
     geno.gp2$magnitude<-as.numeric(as.vector(geno.gp2$magnitude))
     for (i in 1:length(cluster.list)) {
       if  (dim(geno.gp2)[1]==0) {col2<-"white"}
       if (genoclass.nona[2]==cluster.list[i]) {col2<-color.list[i]}
     }
     points(geno.gp2$Angle,geno.gp2$magnitude,col=col2,pch=16,cex=0.8)
     legend(-0.1,-1, legend=c("                             "), bty="o", bg="white",box.col="white", xpd=TRUE, horiz=TRUE)
     legend(-0.1,-1 , legend = c(genoclass.nona[1],genoclass.nona[2],"NA"), col = c(col1, col2, "grey") , pch = 19,cex = 0.9, bty="n", xpd=TRUE, horiz=TRUE)
     if (length(genoclass.nona)>2) { # case with 3 nona clusters
      geno.gp3<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[3]),]
      geno.gp3$Angle<-as.numeric(as.vector(geno.gp3$Angle))
      geno.gp3$magnitude<-as.numeric(as.vector(geno.gp3$magnitude))
      for (i in 1:length(cluster.list)) {
       if  (dim(geno.gp3)[1]==0) {col3<-"white"}
       if (genoclass.nona[3]==cluster.list[i]) {col3<-color.list[i]}
      }
      points(geno.gp3$Angle,geno.gp3$magnitude,col=col3,pch=16,cex=0.8)
      legend(-0.1,-1, legend=c("                            "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)
      legend(-0.1,-1, legend = c(genoclass.nona[1],genoclass.nona[2],genoclass.nona[3],"NA"), col = c(col1, col2, col3,"grey") , pch = 19,cex = 0.9,bty="n", xpd=TRUE, horiz=TRUE)
     } # end condition 3 clusters
    } # end >1 group with no NA genotypes
   }  # end >= 1 cluster which is not NA
  } # end >= 1 cluster  which may be NA
 }  # end condition > 2 called datapoints
dev.off()     # to close plot device for R script in command line
}          # end loop across SNPs

#####################################################################################################################
### Part 3 for hierarchical clustering method providing new clusters plots & genotype calling
# reminder, call.filter1 structure already created whether or not filter 1 applied
# Hclust method will be applied if # hclust.condition<-"yes"  # alternative is "No" and then lines below are not used
# if Hclust method is being asked & filter 1 is skipped (set to 0), it is reset by default to 6 because filter 1 is required to perform the hierarchical clustering
# so to perform Hclust with a chosen filter value, filter 1 needs to be set up to that value

if (hclust.condition=="yes") {
  ## Initializing dataset to be created with Hclust groups info
   snpres.hclust.allchosen.snps<-as.data.frame(matrix(NA,nrow=0,ncol=7,dimnames=list(NULL,c(colnames(snpres.f1),"tree.coord"))))


for (snp in SNPlist){  # loop across SNPs with NA in ori.genotypeId still possible based on Angle & magnitude
    tiff(filename=paste(HclustPlots,"/",snp,".tiff",sep=""),width=690,height=690)     

       snpres.call.1<-as.data.frame(snpres.f1[which(snpres.f1[,"assayId"]==snp),])
       angle.only<-subset(snpres.call.1,select=c(Angle))       # df structure kept if subset data
       rownames(angle.only)<- snpres.call.1$sampleId.new       
       d.angle<-dist(angle.only, method = "euclidean", diag = TRUE, upper = TRUE)   # creates a dist class object
       if (length(snpres.call.1$sampleId.new)<2) {   # if only 1 sample, Clustering not done
         plot(0,0,col="white",pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("less than 2 called datapoints for"," ",snp,sep=""))
       }

       if (length(snpres.call.1$sampleId.new)>=2) {  # if at least 2 called datapoints
           # Clustering with "ward" algorithm
           hclust.angle <- hclust(d.angle, method = "ward.D2")
        # both "ward" and "ward.D" in previous versions are not the correct Ward algo/method , needs "ward.D2" available for R ver 3.1.1  at least

    ###############################################################################################################
    ## Issue of determining how many groups (genoclass) to use for Hclust plot (ie if we have only 2 (or 1) groups)
    ## --> decision to use a priori information < GenotypeId to decide in how many groups to cut the tree

      kgroup<-length(levels(as.factor(as.vector(snpres.call.1$ori.genotypeId))))
      # selecting only genotypes for the particular SNP ifnot all genotypes across snps are inherited
      # NB: if NA data above go into the call and should have been into one of the 3 groups, one more cluster is added
      # which if not possible for biallelic markers--> in this case, Hclust better, but kgroup needs to be 3 only
      #--> add condition below before tree.coord..., also covers case of hand made error by user on raw data
   if (kgroup==4) {kgroup<-3}
   if (kgroup<=3) {
    tree.coord <-cutree(hclust.angle, k=kgroup)
    #  Gets attributions for each individual/sample to the kgroup clusters (genotypes) based on cutting the tree
    #  --> a maximum of 3 clusters will be defined OK, if some individuals are to far away from group, they are retrieved
    # with neighbour distance criterion + warning: each group attribution should have a correspondance with a sampleId new
    tree.coord.df<-cbind(as.data.frame(tree.coord),rownames(as.data.frame(tree.coord)))  #tree.coord.df[1:10,]
    colnames(tree.coord.df)[2]<-"sampleId.new"
    snpres.hclust<-merge(snpres.call.1,tree.coord.df,by.x="sampleId.new",by.y="sampleId.new",all.x=T,all.y=T, sort=T)
    snpres.hclust$Angle<-as.numeric(as.vector(snpres.hclust$Angle))
    snpres.hclust$magnitude<-as.numeric(as.vector(snpres.hclust$magnitude))
    ## defining groups --> if less than 3 groups, gp2 or gp3 will be NULL
    gp1.noNA <-snpres.hclust[which(snpres.hclust$tree.coord=="1" & snpres.hclust$ori.genotypeId!="NA"),]
    gp2.noNA <-snpres.hclust[which(snpres.hclust$tree.coord=="2" & snpres.hclust$ori.genotypeId!="NA"),]
    gp3.noNA <-snpres.hclust[which(snpres.hclust$tree.coord=="3" & snpres.hclust$ori.genotypeId!="NA"),]
    gp1.NA<-snpres.hclust[which(snpres.hclust$tree.coord=="1" & snpres.hclust$ori.genotypeId=="NA"),]
    gp2.NA<-snpres.hclust[which(snpres.hclust$tree.coord=="2" & snpres.hclust$ori.genotypeId=="NA"),]
    gp3.NA<-snpres.hclust[which(snpres.hclust$tree.coord=="3" & snpres.hclust$ori.genotypeId=="NA"),]

    gp1<-snpres.hclust[which(snpres.hclust$tree.coord=="1"),]
    gp2<-snpres.hclust[which(snpres.hclust$tree.coord=="2"),]
    gp3<-snpres.hclust[which(snpres.hclust$tree.coord=="3"),]

## PLots of hierarchical clustering
## defining colors for clusters
cluster.list<-c("A","G","C","T","AG","GA","AC","CA","AT","TA","GC","CG","GT","TG","CT","TC")
color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"black","black","black","black","black","black","black","black","black","black","black","black") # new color list where all heteroz are black   01/04/15

  # reinitialisation geno1/2/3
    geno1<-"  "
      geno2<-"  "
       geno3<-"  "  

for (i in 1:length(cluster.list)) {
    if  (dim(gp1)[1]==0) {col1<-"white"}         
    if  (dim(gp1)[1]!=0) {
      gp1.dfcount<-as.data.frame(table(gp1$ori.genotypeId,gp1$tree.coord))
      most.freqt.Geno1<-gp1.dfcount[which(gp1.dfcount$Freq==max(gp1.dfcount$Freq)),"Var1"]
      rep.type<- length(snpres.hclust[which(snpres.hclust$tree.coord=="1"),which(colnames(snpres.hclust)=="tree.coord")])
      geno1<- as.vector(most.freqt.Geno1[1])
      snpres.hclust[which(snpres.hclust$tree.coord=="1"), which(colnames(snpres.hclust)=="tree.coord")]<- rep(geno1,rep.type)
      if (as.vector(most.freqt.Geno1[1])==cluster.list[i])  {col1<-color.list[i]}
    }
    if  (dim(gp2)[1]==0) {col2<-"white"}
    if (dim(gp2)[1]!=0) { # start loop gp2 # 0
      gp2.dfcount<-as.data.frame(table(gp2$ori.genotypeId,gp2$tree.coord))
      most.freqt.Geno2<-gp2.dfcount[which(gp2.dfcount$Freq==max(gp2.dfcount$Freq)),"Var1"]
      geno2<- as.vector(most.freqt.Geno2[1])
      if  (geno2!="NA") {
        # modif to replace snpres.hclust$tree.coord for group 2 with most frequent allele type of that group i.e as.vector(most.freqt.Geno2[1]) 01/04/15
        # so extract vector within that df corresponding to all "2" value & replace them by the most frequent allele type
        rep.type<- length(snpres.hclust[which(snpres.hclust$tree.coord=="2"),which(colnames(snpres.hclust)=="tree.coord")])
        snpres.hclust[which(snpres.hclust$tree.coord=="2"), which(colnames(snpres.hclust)=="tree.coord")]<- rep(geno2,rep.type)
        if (as.vector(most.freqt.Geno2[1])==cluster.list[i])  {col2<-color.list[i]}
      } # end loop gp2 no NA only
      if (geno2=="NA") {  # start loop if some points have been attributed to Ward clusters that were previously NA
        geno2<-"WC"  # new Ward Cluster identified
        rep.type<- length(snpres.hclust[which(snpres.hclust$tree.coord=="2"),which(colnames(snpres.hclust)=="tree.coord")])
        snpres.hclust[which(snpres.hclust$tree.coord=="2"), which(colnames(snpres.hclust)=="tree.coord")]<- rep(geno2,rep.type)
        col2<-"purple"
      } # end loop new Ward clusters
    } # end loop gp2 # 0
            
    if (dim(gp3)[1]==0) {col3<-"white"}
    if (dim(gp3)[1]!=0) {
      gp3.dfcount<-as.data.frame(table(gp3$ori.genotypeId,gp3$tree.coord))
      most.freqt.Geno3<-gp3.dfcount[which(gp3.dfcount$Freq==max(gp3.dfcount$Freq)),"Var1"]
      geno3<- as.vector(most.freqt.Geno3[1])
      if  (geno3!="NA") {
      # modif to replace snpres.hclust$tree.coord for group 3 with most frequent allele type of that group i.e as.vector(most.freqt.Geno3[1]) 01/04/15
      # so extract vector within that df corresponding to all "3" value & replace them by the most frequent allele type
        rep.type<- length(snpres.hclust[which(snpres.hclust$tree.coord=="3"),which(colnames(snpres.hclust)=="tree.coord")])
        snpres.hclust[which(snpres.hclust$tree.coord=="3"), which(colnames(snpres.hclust)=="tree.coord")]<- rep(geno3,rep.type)
        if (as.vector(most.freqt.Geno3[1])==cluster.list[i])  {col3<-color.list[i]} # so if it's an hetero, it will be black , ifnot see set colours
      } # end loop gp3 no NA only
      if (geno3=="NA") {  # start loop if some points have been attributed to Ward clusters that were previously NA
         geno3<-"WC"  # new Ward Cluster identified
         rep.type<- length(snpres.hclust[which(snpres.hclust$tree.coord=="3"),which(colnames(snpres.hclust)=="tree.coord")])
         snpres.hclust[which(snpres.hclust$tree.coord=="3"), which(colnames(snpres.hclust)=="tree.coord")]<- rep(geno3,rep.type)
         col3<-"purple"
      } # end loop new Ward clusters
    } # end loop gp3 # 0
}  # end loop across all possible types of geno classes

## defining scale 
 if (samescale=="no") {
   mxlim<-max(as.numeric(as.vector(snpres.hclust$Angle)))+0.3
   mylim<-max(as.numeric(as.vector(snpres.hclust$magnitude))+3)
 }
 if (samescale=="yes") {
 mxlim<-max(as.numeric(as.vector(snpres.df$Angle))+0.1)
 mylim<-max(as.numeric(as.vector(snpres.df$magnitude))+3)
 }
 plot(gp1.noNA$Angle,gp1.noNA$magnitude,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(-0.1*mylim,mylim),xlim=c(-0.15,mxlim))
 title(main="Hierarchical Clustering (Ward algorithm)",sub=paste("SNP= ",snp))
 points(gp2.noNA$Angle,gp2.noNA$magnitude,col=col2,pch=16,cex=0.8)
 points(gp3.noNA$Angle,gp3.noNA$magnitude,col=col3,pch=16,cex=0.8)
 points(gp1.NA$Angle,gp1.NA$magnitude,col=col1,pch=10,cex=0.8)
 points(gp2.NA$Angle,gp2.NA$magnitude,col=col2,pch=10,cex=0.8)
 points(gp3.NA$Angle,gp3.NA$magnitude,col=col3,pch=10,cex=0.8)
 legend(-0.1,-1, legend=c("                            "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)
 legend(-0.1,-1, legend = c(geno1,geno2,geno3), col = c(col1, col2, col3) , pch = 19,cex = 0.9,bty="n", xpd=TRUE,horiz=TRUE)
 snpres.hclust.allchosen.snps<-rbind(snpres.hclust.allchosen.snps,snpres.hclust)
  } # end condition for excluding more than 3 groups
 }  # end condition at least 2 called datapoints
 dev.off() 
 } # end loop across snps
} # end hclust.condition=="yes"

 if (hclust.condition=="yes") { # saving of df with hclust groups info
 colnames(snpres.hclust.allchosen.snps)[7]<-"HClust.groups"
 write.table(snpres.hclust.allchosen.snps,file=paste(output.snpcall,"/","res.with.hclust.txt",sep=""),append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)
 } # end hclust.condition == yes for df saving







