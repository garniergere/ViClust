##############################################################################
###### How to cite ViClust:
# Bouteiller XP, Barraquand F, Garnier-G�r� P, Harmand N, Laizet Y, Raimbault A, Segura R, Lassois L, Monty A, Verdu C, Mariette S & Port� AJ (2018)
# No evidence for genetic differentiation in juvenile traits between Belgian and French populations of the invasive tree Robinia pseudoacacia .
# Plant Ecology and Evolution 151 (1): 5�17
## See also Website https://github.com/garniergere/ViClust/ for detailed description,  examples and R editor or standalone version

#############################################################################
# Standalone R script
# deals with raw xml files from Sequenom SNP genotypes
# allows batch visualization and alternative clustering
# Contact: Pauline Garnier-G�r�, pauline.garnier-gere@inra.fr
# March 2018
#############################################################################

###loading libraries & Initialization
library(XML)
##########################################################################
### To launch R script from a command line

cmd_args <- commandArgs(trailingOnly = TRUE)    

## Defining arguments of command line in Rscript
xmlfile<-cmd_args[1] # argument 1 is the xml file to import
snplist<-cmd_args[2] # argument 2 is the chosen SNP list file to plot, put all if all SNP are being analyzed/plotted
indlist<-cmd_args[3] # argument 3 is the individual list file to plot, put all if all individuals are being analyzed/plotted
colclusterA<-cmd_args[4] # colour cluster base "A"    red
colclusterG<-cmd_args[5] # colour cluster base "G"    green
colclusterC<-cmd_args[6] # colour cluster base "C"    blue
colclusterT<-cmd_args[7] # colour cluster base "T"    orange
samescale<-cmd_args[8]   # "yes" or "no" for same scale all plots or not no by default

SequenomPlots <-cmd_args[9] # folder for sequenom plots
output.snpcall<-cmd_args[10] # filename for cluster call results for selected SNPs and individuals  (easier than folder for Galaxy)
html_file_path <- cmd_args[11] # Galaxy html file   --> for each plot to create the html link to view plot
htmlFolder<- cmd_args[12]    # Galaxy html folder   --> is needed to store html files (links) & 3 # folders for all plots

filter1Plots<-cmd_args[13]   #  folder for plots after filter 1
filter1val<-cmd_args[14]     # Minimal magnitude threshold value based on SNR , excluding data with too low signals, 
                             # will plot only called data, default is 6, alternative no or more stringent value  
res.filter1<-cmd_args[15]    # filename for cluster call results after filter 1, skip if filter1val is "0"

hclust.condition<-cmd_args[16]   # "yes" or "no" for performing Hclust method, "no" by default
HclustPlots<-cmd_args[17]    #  folder for plots after Hclust (Hierarchical Clustering)
res.with.hclust<-cmd_args[18] # filename with additional info from Hierarchical clustering after filter 1 (performed with magnitude equals 6 by default if the filter1val is "0")

##Create directory to store images
dir.create(htmlFolder)
SequenomPlotsPath <- paste(htmlFolder, SequenomPlots ,sep="/")
dir.create(SequenomPlotsPath)   # create all path up to the html Folder
filter1PlotsPath <- paste(htmlFolder, filter1Plots ,sep="/")
dir.create(filter1PlotsPath)
HclustPlotsPath <- paste(htmlFolder, HclustPlots ,sep="/")
dir.create(HclustPlotsPath)


############################################################################################################
## Loading xml output/result file
############################################################################################################

data.xml <- xmlParse(xmlfile)   # code with argument 1 in command line  
#data.xml <- xmlParse("EG_OAKTRACK_W1_2_G0772322.xml") # to activate with Tinn R / test with original typeranalyzer & all-records nodes
# data.xml <- xmlParse("Mariette_El_Khoury_Oeno.xml")  # for test 26/03/2015   + NPO setting working dir from R console C:\15-02-28-papers-later+\15-03-Noemie-Stage+portageGalaxy\15-03-Portage-Galaxy\14-11-test-Rscript-Prog

data.list<-xmlToList(data.xml)
records<-t(as.data.frame(data.list$typeranalyzer$'all-records')) ## OK if put all-records after $ in simple quotes ''
#colnames(records)
# [1] "area"            "areaUncert"      "assayId"         "assayPk"         "calibration"     "callPk"          "callScore"       "description"    
# [9] "entryOperator"   "frequency"       "frequencyUncert" "genotypeId"      "height"          "mass"            "maxShift"        "noiseStdDev"    
#[17] "peakScore"       "peakType"        "rasters"         "resolution"      "sampleId"        "samplePk"        "snr"             "spectraPk"      
#[25] "status"          "use-call"        "use-well"        "wellPosition"    "wellsPk"        

record2<-records[,c("sampleId","assayId","genotypeId","height","area","mass","peakScore","peakType","snr","resolution","description","wellPosition")]  # at this stage,w1.oak.rec2 is a matrix, all records values are character 

##########################################################################
### defining list of SNPs and list anf samples (individuals) as macro-variables needed in loop, some refering to arguments above
##########################################################################

### defining list of selected SNP to plot : SNPlist
if (snplist=="all"){  # snplist is argument 2
 SNPlist<-levels(as.factor(record2[,"assayId"]))
}
if (snplist!="all"){
 SNPlist<-scan(snplist,sep="", strip.white=TRUE,what="")
}

### defining list of selected individuals to plot: INDlist

if (indlist=="all"){
INDlist<-levels(as.factor(record2[,"sampleId"]))
}
if (indlist!="all"){
 INDlist<-scan(indlist,sep="", strip.white=TRUE,what="") 
}

### OK to deal with ID unique later when individuals and snps have been selected from initial data, prog will run quicker, see below
### Extract data for chosen SNPs & individuals

selsnp<-matrix(NA,nrow=0,ncol=dim(record2)[2])
colnames(selsnp)<-colnames(record2) 

for (snp in SNPlist) {
     if (snplist=="all")  {selsnp<-record2}                                              # used
     if (snplist!="all") {
           selsnp.s<-record2[which(record2[,"assayId"]==snp),]    
           selsnp<-rbind(selsnp,selsnp.s)
          }
        } 
        
selsnpf<-matrix(NA,nrow=0,ncol=dim(selsnp)[2])
colnames(selsnpf)<-colnames(selsnp) 

for (samp in INDlist)  {  # see indlist argument
     if (indlist=="all") {selsnpf<-selsnp}
     if (indlist!="all") {
           selsnp.i<-selsnp[which(selsnp[,"sampleId"]==samp),]
           selsnpf<-rbind(selsnpf,selsnp.i)
          }
        } 
selsnp<-selsnpf  # final selection of SNP & individuals
colnames(selsnp)[3]<-"ori.genotypeId"

################################################################################
### export raw output from parsed xml file & selected SNP and individual lists
################################################################################
# replace all "" ori.genotypeId by "NA" for exportation to text file without quote
selsnp.export<-selsnp
for (i in 1:dim(selsnp.export)[1]){
 if (selsnp.export[i,"ori.genotypeId"]=="") {selsnp.export[i,"ori.genotypeId"]<-"NA"} 
 }
write.table(as.data.frame(selsnp.export),file=output.snpcall,append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)
### dealing with individuals which are replicated for computing angle & magnitude
sampleId.new<-matrix(NA,ncol=1,nrow=dim(selsnp)[1])
colnames(sampleId.new)<-"sampleId.new"
for (i in 1:dim(selsnp)[1]) {
sampleId.new[i,1]<-paste(selsnp[i,"sampleId"],".",selsnp[i,"wellPosition"],sep="")
}
selsnp1<-selsnp
selsnp<-cbind(sampleId.new,selsnp1)        
# selsnp[1:10,]

#############################################
## create Matrix with one line per combination SNP*individual*well (sampleId.new, see above) & computing angle value, magnitude value + NA instead of ""
snpres<-matrix(NA,nrow=0,ncol=5)
colnames(snpres)<-c(colnames(selsnp)[c(1,3,4)],"Angle","magnitude") 

# defining complete sample list of non redundant individuals
samplelist<-levels(as.factor(selsnp[,"sampleId.new"]))

for (snp in SNPlist) {
     for (sam in samplelist) {
      samp<-selsnp[which(selsnp[,"assayId"]==snp & selsnp[,"sampleId.new"]==sam),]  
      if (length(samp)!=0) {   # for non-empty snpres files
      highmass.record<-samp[which(samp[,"mass"]==max(samp[,"mass"])),]  
      samp.less.probe<-samp[-which(samp[,"mass"]==min(samp[,"mass"])),]
      lowmass.record<-samp.less.probe[which(samp.less.probe[,"mass"]==min(samp.less.probe[,"mass"])),]
      if (as.numeric(lowmass.record["height"])==0) {
      samp.res<-t(matrix(c(sam,snp,samp[1,"ori.genotypeId"],pi/2,
                sqrt((as.numeric(highmass.record["snr"]))^2 + (as.numeric(lowmass.record["snr"]))^2)),dimnames=list(colnames(snpres))))}
       if (as.numeric(lowmass.record["height"])!=0) {
      samp.res<-t(matrix(c(sam,snp,samp[1,"ori.genotypeId"],atan(as.numeric(highmass.record["height"])/as.numeric(lowmass.record["height"])), # use height for angle
                  sqrt((as.numeric(highmass.record["snr"]))^2 + (as.numeric(lowmass.record["snr"]))^2)),dimnames=list(colnames(snpres)))) # but snr for magnitude
                    }
         if (samp.res[,"ori.genotypeId"]=="") {samp.res[,"ori.genotypeId"]<-"NA"} # --> line works if we keep matrix structure but not the one with is.na    
         snpres<-rbind(snpres,samp.res)  
          }   # end loop for non empty snpres file
        }     # end loop samples (or individuals) with new idividuals ID (no replicate ID)
      } # end SNPlist

############################################################################################################
## Loop 1 for plotting SNPs (& assigning them to clusters) as in the Sequenom sofware after the Cluster Call
############################################################################################################

################################
#Create html file listing images
################################
# setting up the html file path list that includes the html files for Sequenom Plots
write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">', file = html_file_path)
write('<html lang="en" xml:lang="en" xmlns="http://www.w3.org/1999/xhtml"><head><body>', file = html_file_path, append = TRUE)
write("<h1>Sequenom Plots</h1><br/>", file = html_file_path, append = TRUE)                    ####  html_file_path is argument 11
write("<p>Click on a link to view a plot</p><br/>", file = html_file_path, append = TRUE)
write("<h2>Sequenom plots</h2><br/>", file = html_file_path, append = TRUE) 

for (snp in SNPlist){  # loop across SNPs
png(filename=paste(SequenomPlotsPath,"/",snp,".png",sep=""),width=690,height=690)  # all path is needed

       snpres.1<-as.data.frame(snpres[which(snpres[,"assayId"]==snp),])  
        if (length(snpres.1$sampleId.new)<2) {   # if only 1 datapoint
         plot(0,0,col="white",pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("less than 2 called datapoints for"," ",snp,sep=""))
         }
       if (length(snpres.1$sampleId.new)>=2) {  # if at least 2 called datapoints  15/11/14
 
       genoclass<-levels(as.factor(as.vector(snpres.1[,"ori.genotypeId"])))
       if (length(genoclass)==0) {   # if no data
         plot(0,0,col="white",pch=16,cex=1, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("No data for"," ",snp,sep=""))}

      if (length(genoclass)!=0) {   # start loop at least one cluster  (including NA cluster)
      if (length(which(genoclass=="NA"))==0) {genoclass.nona<-genoclass}   
      if (length(which(genoclass=="NA"))!=0)
      {genoclass.nona<-genoclass[-which(genoclass=="NA")]} # remove "NA" levels if they do exist
      
### Plots different group of genotypes with colours as argument in command line
## samescale<-"no"
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
#color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"purple","purple","purple2","purple2","purple4","purple4","cyan","cyan","turquoise2","turquoise2","turquoise4","turquoise4")
color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"black","black","black","black","black","black","black","black","black","black","black","black") # new color list where all heteroz are black   01/04/15

for (i in 1:length(cluster.list)) {
    if  (dim(geno.gp)[1]==0) {col1<-"white"}      
    else if (genoclass.nona[1]==cluster.list[i]) {col1<-color.list[i]}
 }

 if (dim(geno.gp)[1]==0) { # if only NA data
     plot(0,0,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,mylim),xlim=c(0,mxlim),main=paste("Only NoCall data for"," ",snp,sep=""))
      # defining & plotting eventual nocall cluster
     geno.nocall<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
     geno.nocall$Angle<-as.numeric(as.vector(geno.nocall$Angle))   
     geno.nocall$magnitude<-as.numeric(as.vector(geno.nocall$magnitude))
     points(geno.nocall$Angle,geno.nocall$magnitude,col="grey",pch=16,cex=0.8)
     }
     
 if (dim(geno.gp)[1]!=0) { #if at least one class which is not NA (condition could have been done with genoclass.nona
      plot(geno.gp$Angle,geno.gp$magnitude,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(-0.1*mylim,mylim),xlim=c(-0.15,mxlim),main=snp)
      legend(-0.1,-1 , legend = c(genoclass.nona[1], "NA"), col = c(col1,"grey"), pch=19, cex=0.9, bty="n", xpd=TRUE, horiz=TRUE) # modif 30/03/15

      # defining & plotting eventual NA cluster
     geno.nocall<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
     geno.nocall$Angle<-as.numeric(as.vector(geno.nocall$Angle))   
     geno.nocall$magnitude<-as.numeric(as.vector(geno.nocall$magnitude))
     points(geno.nocall$Angle,geno.nocall$magnitude,col="grey",pch=16,cex=0.8)

  if (length(genoclass.nona)>1) { # if more than one class with called genotypes
     geno.gp2<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[2]),]   
     geno.gp2$Angle<-as.numeric(as.vector(geno.gp2$Angle))   
     geno.gp2$magnitude<-as.numeric(as.vector(geno.gp2$magnitude))
     for (i in 1:length(cluster.list)) {
      if  (dim(geno.gp2)[1]==0) {col2<-"white"}    # added 14/10/12
      if (genoclass.nona[2]==cluster.list[i]) {col2<-color.list[i]}
     } 
     points(geno.gp2$Angle,geno.gp2$magnitude,col=col2,pch=16,cex=0.8)
     # 2 lines below modif 30/03/15
     legend(-0.1,-1, legend=c("                             "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)  # draws a white rectangle to cover previous legend
    legend(-0.1,-1 , legend = c(genoclass.nona[1],genoclass.nona[2],"NA"), col = c(col1, col2, "grey") , pch = 19,cex = 0.9, bty="n", xpd=TRUE, horiz=TRUE)
        
     if (length(genoclass.nona)>2) { # case with 3 nona clusters
      geno.gp3<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[3]),]
      geno.gp3$Angle<-as.numeric(as.vector(geno.gp3$Angle))   
      geno.gp3$magnitude<-as.numeric(as.vector(geno.gp3$magnitude))
      for (i in 1:length(cluster.list)) {
        if  (dim(geno.gp3)[1]==0) {col3<-"white"}    # added 14/10/12
        if (genoclass.nona[3]==cluster.list[i]) {col3<-color.list[i]}
      }
      points(geno.gp3$Angle,geno.gp3$magnitude,col=col3,pch=16,cex=0.8)
      # 2 lines below : modif 30/03/15
    legend(-0.1,-1, legend=c("                            "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)
    legend(-0.1,-1, legend = c(genoclass.nona[1],genoclass.nona[2],genoclass.nona[3],"NA"), col = c(col1, col2, col3,"grey") , pch = 19,cex = 0.9,bty="n", xpd=TRUE, horiz=TRUE)

   } # end condition if 3 clusters of genotypes with or without nocall NA
  }  # end loop more than 1 cluster call of genotypes to plot with or without nocall NA
 } # end loop at least 1 called cluster with or without no call points
} # end loop at least 1 cluster (including no call)
write(paste('<a href="', paste(SequenomPlots,"/",snp,".png",sep=""), '">Sequenom Plot for ', snp, '</a><br/>', sep=""), file = html_file_path, append = TRUE)
         # line above for saving the plot accessible in an html file from a list defined by html_file_path, one html file for each snp so in the loop
         # a html file is like an xml file with nodes # modif syntax 29/10/2014
 }  # end condition at least 2 called datapoints  15/11/14
dev.off()   
}          # end loop across SNPs & SEQUENOM clusters

####################################################################################################
#### Loop 2 for filtering based on magnitude and plotting data after excluding filtered points 
####################################################################################################

#### Create a call factor for excluding datapoints with too low magnitude based on SNR (Filter 1)
  call.filter1<-matrix(NA,ncol=1,nrow=dim(snpres)[1])
  colnames(call.filter1)<-"call.filter1"

#filter1val<-6   # argument 14 put arg value in command line directly, default is "6": Minimal magnitude threshold value for peaks to be kept 
## Lines below added for exporting Angle & Magnitude in both cases 15/11/13 : creation snpres2 with col nocall.below6 / cal.above6 + print snpres with angle & magnitude

#head(snpres)
 snpres.df<-as.data.frame(snpres)  #  is.data.frame(snpres.df)
#head(snpres.df)
# is.factor(snpres.df$magnitude)
 snpres.df$magnitude<-as.numeric(as.vector(snpres.df$magnitude))     
#head(snpres.df)

#filter1val<-0  
filter1val<-as.numeric(filter1val)
# output.snpcall<-"output_Files14"

if (filter1val==0) { 
  for (i in 1:dim(snpres.df)[1]) {
  # i<-1  
  if (snpres.df$magnitude[i]<=6) {call.filter1[i,1]<-"nocall.below6" }
  else if (snpres.df$magnitude[i]>6) {call.filter1[i,1]<-"call.above6" }
  }   
  snpres.df2<-cbind(snpres.df,call.filter1)       
  snpres2<-snpres.df2 # to adjust to code below                        
write.table(as.data.frame(snpres2),file=res.filter1,append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE) # for Galaxy
  } # end condition for applying no filter 1 a priori --> file with Angle & Magnitude is still created with call column for magnitude 6

if (filter1val>0) {  # so if filter1val=0, see above
 for (i in 1:dim(snpres.df)[1]) {
  # i<-1  
  if (snpres.df$magnitude[i]<=filter1val) {call.filter1[i,1]<-"nocall" }  
  else if (snpres.df$magnitude[i]>filter1val) {call.filter1[i,1]<-"call" }
  } 
 snpres.df2<-cbind(snpres.df,call.filter1)       
 snpres2<-snpres.df2 
######## export file with ori.genotypeId + information "call" & "nocall"--> will go directly into Galaxy folder
write.table(as.data.frame(snpres2),file=res.filter1,append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE) # res.filter1=argument 15
} # end condition filter1                                            

######################################################################                                                                     
########## Plots if using filter 1 based on SNR ######################
########## should include all grey data corresponding to "nocall" based on "" ori.genotypeId
   
##  creation snpres.f1 in all cases filter1val before loops for Plots
 if (filter1val==0) { #start condition filter 1
    snpres.f1<-subset(snpres2,subset=call.filter1=="call.above6")    # Df with only "call.above6" points based on first filter on magnitude
} # end condition filter 1 

if (filter1val>0) { #start condition filter 1
    snpres.f1<-subset(snpres2,subset=call.filter1=="call")          # Dataframe with only "call" points based on first filter on magnitude
} # end condition filter 1 


#########################################
### Loop on plots after applying filter 1

## Appending new html files (plots) to the list previously defined
write("<h2>Filter 1 plots</h2><br/>", file = html_file_path, append = TRUE)

for (snp in SNPlist){  # loop across SNPs still NA possible in ori.genotypeId
 png(filename=paste(filter1PlotsPath,"/",snp,".png",sep=""),width=690,height=690) 

# tiff(filename=paste(filter1PlotsPath,"/",snp,".tiff",sep=""),width=690,height=690) 
   
       snpres.1<-as.data.frame(snpres.f1[which(snpres.f1[,"assayId"]==snp),])  # change snpres in snpres.f1 for this loop 
            if (length(snpres.1$sampleId.new)<2) {   # if only 1 sample
         plot(0,0,col="white",pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("less than 2 called datapoints for"," ",snp,sep=""))
       #  savePlot(filename=paste("HCLUST.PLOTs.test141015-ward.D","/",snp,".tiff",sep=""), type="tiff")
         }
       if (length(snpres.1$sampleId.new)>=2) {  # if at least 2 called datapoints  15/11/14

       genoclass<-levels(as.factor(as.vector(snpres.1[,"ori.genotypeId"])))

       if (length(genoclass)==0) {   # if no data
         plot(0,0,col="white",pch=16,cex=1, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("No data for"," ",snp,sep=""))}

       if (length(genoclass)!=0) {   # start loop at least one cluster  which may be NA
       if (length(which(genoclass=="NA"))==0) {genoclass.nona<-genoclass}
       if (length(which(genoclass=="NA"))!=0)
       {genoclass.nona<-genoclass[-which(genoclass=="NA")]} # remove "NA" levels if they do exist  in nona group/cluster

### Plots different group of genotypes with colours as argument in command line
## samescale<-"no"
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
#color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"purple","purple","purple2","purple2","purple4","purple4","cyan","cyan","turquoise2","turquoise2","turquoise4","turquoise4")
color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"black","black","black","black","black","black","black","black","black","black","black","black") # new color list where all heteroz are black   01/04/15

for (i in 1:length(cluster.list)) {
    if  (dim(geno.gp)[1]==0) {col1<-"white"}      
    else if (genoclass.nona[1]==cluster.list[i]) {col1<-color.list[i]}
   }

   if (dim(geno.gp)[1]==0) { # if only NA data
     plot(0,0,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,mylim),xlim=c(0,mxlim),main=paste("Only NoCall data for"," ",snp,sep=""))
      # defining & plotting eventual nocall cluster
     geno.nocall<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
     geno.nocall$Angle<-as.numeric(as.vector(geno.nocall$Angle))
     geno.nocall$magnitude<-as.numeric(as.vector(geno.nocall$magnitude))
     points(geno.nocall$Angle,geno.nocall$magnitude,col="grey",pch=16,cex=0.8)
     }
  if (dim(geno.gp)[1]!=0) { #if at least one class which is not NA (condition could have been done with genoclass.nona
      plot(geno.gp$Angle,geno.gp$magnitude,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(-0.1*mylim,mylim),xlim=c(-0.15,mxlim),main=snp)
      legend(-0.1,-1 , legend = c(genoclass.nona[1], "NA"), col = c(col1,"grey"), pch=19, cex=0.9, bty="n", xpd=TRUE, horiz=TRUE)
 
      # defining & plotting eventual NA cluster
     geno.NA<-snpres.1[which(snpres.1$ori.genotypeId=="NA"),]
     geno.NA$Angle<-as.numeric(as.vector(geno.NA$Angle))
     geno.NA$magnitude<-as.numeric(as.vector(geno.NA$magnitude))
     points(geno.NA$Angle,geno.NA$magnitude,col="grey",pch=16,cex=0.8)

  if (length(genoclass.nona)>1) { # if more than one class with no NA  genotypes
     geno.gp2<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[2]),]
     geno.gp2$Angle<-as.numeric(as.vector(geno.gp2$Angle))   
     geno.gp2$magnitude<-as.numeric(as.vector(geno.gp2$magnitude))

     for (i in 1:length(cluster.list)) {
       if  (dim(geno.gp2)[1]==0) {col2<-"white"}    # added 14/10/12
       if (genoclass.nona[2]==cluster.list[i]) {col2<-color.list[i]}
     }
     points(geno.gp2$Angle,geno.gp2$magnitude,col=col2,pch=16,cex=0.8)
     # 2 lines below : modif 30/03/15
    legend(-0.1,-1, legend=c("                            "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)  # draws a white rectangle to cover previous legend
    legend(-0.1,-1 , legend = c(genoclass.nona[1],genoclass.nona[2],"NA"), col = c(col1, col2, "grey") , pch = 19,cex = 0.9, bty="n", xpd=TRUE, horiz=TRUE)
         
     if (length(genoclass.nona)>2) { # case with 3 nona clusters
     geno.gp3<-snpres.1[which(snpres.1$ori.genotypeId==genoclass.nona[3]),]
     geno.gp3$Angle<-as.numeric(as.vector(geno.gp3$Angle))   
     geno.gp3$magnitude<-as.numeric(as.vector(geno.gp3$magnitude))
     for (i in 1:length(cluster.list)) {
      if  (dim(geno.gp3)[1]==0) {col3<-"white"}    # added 14/10/12
      if (genoclass.nona[3]==cluster.list[i]) {col3<-color.list[i]}
      }
      points(geno.gp3$Angle,geno.gp3$magnitude,col=col3,pch=16,cex=0.8)
   legend(-0.1,-1, legend=c("                             "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)
   legend(-0.1,-1, legend = c(genoclass.nona[1],genoclass.nona[2],genoclass.nona[3],"NA"), col = c(col1, col2, col3,"grey") , pch = 19,cex = 0.9,bty="n", xpd=TRUE, horiz=TRUE)

   
    } # end condition 3 clusters
   } # end more than one class/group with no NA genotypes
  }  # end loop at least 1 cluster which is not NA
 } # end loop at least 1 cluster  which may be NA

 write(paste('<a href="', paste(filter1Plots,"/",snp,".png",sep=""), '">Filter 1 Plot for ', snp, '</a><br/>', sep=""), file = html_file_path, append = TRUE)
 }  # end condition at least 2 called datapoints  15/11/14
dev.off()  
  }          # end loop across SNPs

#######################################################################################
### Loop 3 for new clusters plots & genotype calling based on hierarchical clustering   
# reminder, call.filter1 structure already created whether or not filter 1 applied

# Hclust method will be applied if # hclust.condition<-"yes"  # alternative is "No" and then lines below are not activated
# if Hclust method asked & filter 1 skipped (set to 0), it is reset by default to 6 because filter 1 is required to perform the hierarchical clustering 
# lines below applied only if hclust.condition is "yes"       
# hclust.condition<-"no"  # default

if (hclust.condition=="yes") {   # lines applied below if hclust wanted  
  ## Initializing final dataset with Hclust groups info
  snpres.hclust.allchosen.snps<-as.data.frame(matrix(NA,nrow=0,ncol=7,dimnames=list(NULL,c(colnames(snpres.f1),"tree.coord")))) 

## Appending new html files (Hclust plots) to the list previously defined
write("<h2>Plots based on hierarchical clustering</h2><br/>", file = html_file_path, append = TRUE)
 
 for (snp in SNPlist){  # loop across SNPs  with NA in ori.genotypeId still possible based on Angle & magnitude

png(filename=paste(HclustPlotsPath,"/",snp,".png",sep=""),width=690,height=690) 
       snpres.call.1<-as.data.frame(snpres.f1[which(snpres.f1[,"assayId"]==snp),])  # snpres.call.1 is for the one selected snp
     # Extracting angle values to create distance matrix d.angle 
       angle.only<-subset(snpres.call.1,select=c(Angle))                
       rownames(angle.only)<- snpres.call.1$sampleId.new       
       d.angle<-dist(angle.only, method = "euclidean", diag = TRUE, upper = TRUE)   # creates a dist class object 
     # Clustering with "ward" algorithm  

     if (length(snpres.call.1$sampleId.new)<2) {   # if only 1 sample, Clustering can't be done condition added 15/10/14
         plot(0,0,col="white",pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(0,10),xlim=c(0,10),main=paste("less than 2 called datapoints for"," ",snp,sep=""))
         write(paste('<a href="', paste(HclustPlots,"/",snp,".tiff",sep=""), '">Hclust Plot for ', snp, '</a><br/>', sep=""), file = html_file_path, append = TRUE)
     }

     if (length(snpres.call.1$sampleId.new)>=2) {  # if at least 2 called datapoints
       hclust.angle <- hclust(d.angle, method = "ward.D2") # change "ward" in "ward.D" for version R 3.1.1, but pb both "ward" and
                                                           # "ward.D" are not the correct ward algo/method , needs "ward.D2" available for R ver 3.1.1  at least
     # Issue of max nb of groups --> decision to keep 3 max & as many groups as for the typer 4 mix of Gaussian method
       kgroup<-length(levels(as.factor(as.vector(snpres.call.1$ori.genotypeId))))
     # Note 14/10/12: if NA data above which go into the call and should have been into one of the 3 groups, it is adding one class ie 4 in total
     # which if not possible, so in this case, Hclust will propose something better, but kgroup needs to be 3 only --> add a condition before tree.coord
    if (kgroup==4) {kgroup<-3} # added condition 14/10/12 & below too
    if (kgroup<=3) {  # condition for excluding "silly" or "user call mistakes" cases 
       tree.coord <-cutree(hclust.angle, k=kgroup)   #  attribute each individual/sample to the kgroup clusters based on tree cutting
     # Merging snpres & Hclust group info
       tree.coord.df<-cbind(as.data.frame(tree.coord),rownames(as.data.frame(tree.coord)))  #tree.coord.df[1:10,]
       colnames(tree.coord.df)[2]<-"sampleId.new"                               
       snpres.hclust<-merge(snpres.call.1,tree.coord.df,by.x="sampleId.new",by.y="sampleId.new",all.x=T,all.y=T, sort=T)  
       snpres.hclust$Angle<-as.numeric(as.vector(snpres.hclust$Angle))  # ifnot plot of factor give boxplot
       snpres.hclust$magnitude<-as.numeric(as.vector(snpres.hclust$magnitude))

     # defining groups --> if less than 3 groups, gp2 or gp3 will be NULL

gp1.noNA <-snpres.hclust[which(snpres.hclust$tree.coord=="1" & snpres.hclust$ori.genotypeId!="NA"),]
gp2.noNA <-snpres.hclust[which(snpres.hclust$tree.coord=="2" & snpres.hclust$ori.genotypeId!="NA"),]
gp3.noNA <-snpres.hclust[which(snpres.hclust$tree.coord=="3" & snpres.hclust$ori.genotypeId!="NA"),]
gp1.NA<-snpres.hclust[which(snpres.hclust$tree.coord=="1" & snpres.hclust$ori.genotypeId=="NA"),]
gp2.NA<-snpres.hclust[which(snpres.hclust$tree.coord=="2" & snpres.hclust$ori.genotypeId=="NA"),]
gp3.NA<-snpres.hclust[which(snpres.hclust$tree.coord=="3" & snpres.hclust$ori.genotypeId=="NA"),]

gp1<-snpres.hclust[which(snpres.hclust$tree.coord=="1"),]
gp2<-snpres.hclust[which(snpres.hclust$tree.coord=="2"),]
gp3<-snpres.hclust[which(snpres.hclust$tree.coord=="3"),]

##### PLots of hierarchical clustering

## defining colors for new Hclust groups
cluster.list<-c("A","G","C","T","AG","GA","AC","CA","AT","TA","GC","CG","GT","TG","CT","TC")
# colclusterA,colclusterG,colclusterC,colclusterT below have been defined in loop sequenom plots above for windows R prog
#color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"purple","purple","purple2","purple2","purple4","purple4","cyan","cyan","turquoise2","turquoise2","turquoise4","turquoise4")
color.list<-c(colclusterA,colclusterG,colclusterC,colclusterT,"black","black","black","black","black","black","black","black","black","black","black","black") # new color list where all heteroz are black   01/04/15

  # reinitialisation geno1/2/3
    geno1<-"  "
      geno2<-"  "
       geno3<-"  "  

for (i in 1:length(cluster.list)) {   # across all possible types of genotypes with no NA in ori.genotypeId for attributing color sets
   # i<-15    
    if  (dim(gp1)[1]==0) {col1<-"white"}         
    if  (dim(gp1)[1]!=0) {
    #    geno1<-"  "
        gp1.dfcount<-as.data.frame(table(gp1$ori.genotypeId,gp1$tree.coord)) #   colnames(gp1)
        most.freqt.Geno1<-gp1.dfcount[which(gp1.dfcount$Freq==max(gp1.dfcount$Freq)),"Var1"]
        
        # modif to replace snpres.hclust$tree.coord for group 1 with most frequent allele type of that group i.e as.vector(most.freqt.Geno1[1]) 01/04/15
        # so extract vector within that df corresponding to all "1" value & replace them by the most frequent allele type
          rep.type<- length(snpres.hclust[which(snpres.hclust$tree.coord=="1"),which(colnames(snpres.hclust)=="tree.coord")])
          geno1<- as.vector(most.freqt.Geno1[1])     
          snpres.hclust[which(snpres.hclust$tree.coord=="1"), which(colnames(snpres.hclust)=="tree.coord")]<- rep(geno1,rep.type)
          
          if (as.vector(most.freqt.Geno1[1])==cluster.list[i])  {col1<-color.list[i]}
      }
    if  (dim(gp2)[1]==0) {col2<-"white"}
    if(dim(gp2)[1]!=0) { # start loop gp2 # 0 
     #   geno2<-"  "
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
             }# end loop new Ward clusters
          } # end loop gp2 # 0 
            
   if  (dim(gp3)[1]==0) {col3<-"white"}
   if (dim(gp3)[1]!=0) { # start loop gp3 # 0 
      #  geno3<-"    "
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
# samescale<-"no"        # "no" by default
 if (samescale=="no") {   
mxlim<-max(as.numeric(as.vector(snpres.hclust$Angle))+0.3)
mylim<-max(as.numeric(as.vector(snpres.hclust$magnitude))+3)}
 if (samescale=="yes") {
 mxlim<-max(as.numeric(as.vector(snpres.df$Angle))+0.1)
 mylim<-max(as.numeric(as.vector(snpres.df$magnitude))+3) }  

# change below all gp1,2,3 in 6 groups separating noNA & NA in ori.genotypeID
plot(gp1.noNA$Angle,gp1.noNA$magnitude,col=col1,pch=16,cex=0.8, xlab="Angle",ylab="Magnitude",ylim=c(-0.1*mylim,mylim),xlim=c(-0.15,mxlim))
title(main="Hierarchical Clustering (Ward algorithm)",sub=paste("SNP= ",snp))
points(gp2.noNA$Angle,gp2.noNA$magnitude,col=col2,pch=16,cex=0.8)
points(gp3.noNA$Angle,gp3.noNA$magnitude,col=col3,pch=16,cex=0.8)
points(gp1.NA$Angle,gp1.NA$magnitude,col=col1,pch=10,cex=0.8)
points(gp2.NA$Angle,gp2.NA$magnitude,col=col2,pch=10,cex=0.8)
points(gp3.NA$Angle,gp3.NA$magnitude,col=col3,pch=10,cex=0.8)
  
   legend(-0.1,-1, legend=c("                         "), bty="o", bg="white", box.col="white", xpd=TRUE, horiz=TRUE)
   legend(-0.1,-1, legend = c(geno1,geno2,geno3), col = c(col1, col2, col3) , pch = 19,cex = 0.9,bty="n", xpd=TRUE,horiz=TRUE)

 write(paste('<a href="', paste(HclustPlots,"/",snp,".png",sep=""), '">Hclust Plot for ', snp, '</a><br/>', sep=""), file = html_file_path, append = TRUE)
 ## Appending df with hierarchical clustering group info
 snpres.hclust.allchosen.snps<-rbind(snpres.hclust.allchosen.snps,snpres.hclust)
   
   } # end condition for excluding more than 3 groups
  }  # end condition at least 2 called datapoints
 dev.off() 
 } # end loop across snps
} # end hclust.condition=="yes"

if (hclust.condition=="yes") { # saving of df with hclust groups information if hclust condition is yes and so lines above have been applied 
  colnames(snpres.hclust.allchosen.snps)[7]<-"HClust.groups"  # head(snpres.hclust.allchosen.snps)
  write.table(snpres.hclust.allchosen.snps,file=res.with.hclust,append=FALSE,quote=FALSE,sep="\t",col.names=TRUE, row.names=FALSE)    # res.with.clust is arg 17
} # end hclust.condition == yes for df saving

## closing html list
write("</body></head></html>", file = html_file_path, append = TRUE)   


