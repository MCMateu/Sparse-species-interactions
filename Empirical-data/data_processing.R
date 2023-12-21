'
Program to compute the covariance matrix for every biome
and for a list of different (abundance) occupancies
'
#remove all variables previously saved in the R environment
rm(list=ls())

##################################
#       LIBRARIES                #
##################################
library(tidyverse)

###################################
#           LOAD DATA             #  
###################################

#Input Folder
InputFolder<-"/home/jose/PCA"
#Output Folder
OutputFolder<-"/home/jose/PCA"


#Load data, data is in long-format
load(paste0(InputFolder,"/dataestimate.Rdata"))

#First Filter
nreads_cutoff <- 1^1
count_cutoff <- 0

remove_runs = c("ERR1104477", "ERR1101508", "SRR2240575") ## duplicate runs for the same sample
DF<- proj %>% filter(  nreads > nreads_cutoff, count > count_cutoff ) %>% filter( ! run_id %in% remove_runs ) %>% 
  group_by( project_id, classification, otu_id ) %>%  mutate(p=count/nreads,o=n()) %>% ungroup() %>% 
  mutate(idall = paste(project_id, classification))

unique(DF$idall)

#Environments to be maintained
names<-c("MGYS00002437  seawater","SRP056641 GUT","SRP052295 Environmental Terrestrial Soil","SRP056641 ORAL",
         "ERP017997 Glacier","ERP012927 River","ERP012927 Lake","ERP009143 activatedsludge" ,"ERP015450 GUT")

DF <- DF %>% filter(idall %in% names)

#Compute occupancy
DF <- DF %>% group_by( project_id, idall) %>% mutate(o=o/n_distinct(run_id))

###############################################################
#Replace the names of the environments for the sake of clarity#
###############################################################

new_names<-c("Sludge", "Lake", "River", "Gut2","Glacier", "Seawater", "Soil", "Gut1", "Oral1")
old_names<- c("ERP009143 activatedsludge","ERP012927 Lake","ERP012927 River","SRP056641 GUT",
              "ERP017997 Glacier","MGYS00002437  seawater",
              "SRP052295 Environmental Terrestrial Soil","ERP015450 GUT","SRP056641 ORAL")

shortnames <- data.frame(new_names,old_names)

#Total number of biomes
N<-length(unique(DF$idall))

#Replace names
for( i in 1:N) {
  
  DF$idall[DF$idall==old_names[i]]<-new_names[i]
  
}  

#List of biomes
biomeList<-unique(DF$idall)

####################################################
#               LOOP OVER OCCUPANCIES              #
####################################################

o_List<-c(0.5,0.95,0.99)

for(o_cutoff in o_List){

    #####################################################
    #               LOOP OVER ALL BIOMES                #
    #####################################################
    for(biome in biomeList){
      
      #Print biome name in the console
      print(biome)
      
      #save single biome data
      biomeDF <- DF %>% filter(idall==biome) %>% ungroup()
      
      #filter by occupancy
      biomeDF <- biomeDF %>% filter(o>o_cutoff)
      
      #Convert wide format
      biomeDF <- spread(data=biomeDF[,c("otu_id","p","sample_id")],key=sample_id,value=p)
      
      #Replace NA by 0
      biomeDF <-replace(biomeDF , is.na(biomeDF ), 0)
      #remove ASV_ID column 
      biomeDF <-subset(biomeDF ,select=-c(otu_id))
    
      #total number of OTUs (species)
      S= nrow(biomeDF)
      
      '
      compute the covariance matrix
      -cov() function, rows: observations, cols: variables
      '
      
      #transpose matrix for correct input to cov() function
      biomeDF<-t(biomeDF)
      
      CovMat<-cov(biomeDF)
      
      PearsMat<-cor(biomeDF, method = "pearson")
      
      '
      Save CovMat in a Folder
      1-Create Abundance Covariance Matriz Folder
      2-Create Biome Folder
      3-Create Occupancy Folder
      4-Save data
      '
      
      #1-Create Abundance Covariance Matriz Folder
      CovFolder<-paste0(OutputFolder,"/AbundanceCov")
      
      if (dir.exists(CovFolder)==FALSE){ dir.create(CovFolder)} else {
        print("Cov Dir already exists!")
      }
      
      #2-Create Biome Folder
      biomeFolder<-paste0(CovFolder,"/",biome)
      if (dir.exists( biomeFolder)==FALSE){ dir.create( biomeFolder)} else {
          print("Cov Dir already exists!")
      }
      
      #3-Create Occupancy Folder
      o_cutoff_Folder<-paste0( biomeFolder,"/o_",o_cutoff)
      
      if (dir.exists(o_cutoff_Folder)==FALSE){ dir.create(o_cutoff_Folder)} else {
        print("Cov Dir already exists!")
      }
      
      write.csv(CovMat,file=paste0(o_cutoff_Folder,"/AbdCovMat.csv"))
      
      '
      Save pearson correlation matrix 
      1-Create Pearson Correlation Abundance Matriz Folder
      2-Create Biome Folder
      3-Create Occupancy Folder
      4-Save data
      '
      
      #1-Create Pearson Correlation Abundance Matriz Folder
      PearsFolder<-paste0(OutputFolder,"/PearsMat")
      
      if (dir.exists(PearsFolder)==FALSE){ dir.create(PearsFolder)} else {
        print("Pears Dir already exists!")
      }
      
      #2-Create Biome Folder
      biomeFolder<-paste0(PearsFolder,"/",biome)
      if (dir.exists( biomeFolder)==FALSE){ dir.create( biomeFolder)} else {
        print("Biome Dir already exists!")
      }
      
      #3-Create Occupancy Folder
      o_cutoff_Folder<-paste0( biomeFolder,"/o_",o_cutoff)
      
      if (dir.exists(o_cutoff_Folder)==FALSE){ dir.create(o_cutoff_Folder)} else {
        print("Occupancy Dir already exists!")
      }
      
      write.csv(PearsMat,file=paste0(o_cutoff_Folder,"/PearsMat.csv"))
      
      '
      Save wide format abudance data 
      1-Create wide Abundance Matriz Folder
      2-Create Biome Folder
      3-Create Occupancy Folder
      4-Save data
      '
      
      #1-Create wide Abundance Matriz Folder
      AbdFolder<-paste0(OutputFolder,"/AbundanceWide")
      
      if (dir.exists(AbdFolder)==FALSE){ dir.create(AbdFolder)} else {
        print("Abundance Dir already exists!")
      }
      
      #2-Create Biome Folder
      biomeFolder<-paste0(AbdFolder,"/",biome)
      if (dir.exists( biomeFolder)==FALSE){ dir.create( biomeFolder)} else {
        print("Biome Dir already exists!")
      }
      
      #3-Create Occupancy Folder
      o_cutoff_Folder<-paste0( biomeFolder,"/o_",o_cutoff)
      
      if (dir.exists(o_cutoff_Folder)==FALSE){ dir.create(o_cutoff_Folder)} else {
        print("Occupancy Dir already exists!")
      }
      
      write.csv(biomeDF,file=paste0(o_cutoff_Folder,"/AbdWide.csv"))
     
    }#END OF THE LOOP OVER BIOMES

}#END OF THE LOOP OVER OCCUPANCIES

