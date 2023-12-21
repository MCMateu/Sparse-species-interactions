rm(list=ls())
##################################
#       LIBRARIES                #
##################################
library(tidyverse)
library(ggplot2)
library(cowplot)
library(Matrix)
library(fitdistrplus)
library(RColorBrewer)
library(scales)
library(pracma)
##################################
#       LIBRARIES                #
##################################

##################################
#       FUNCTIONS                #
##################################
source("/home/jose/MEGA/Figures/Functions.R")
##################################
#       FUNCTIONS                #
##################################

#Main Folder where the simulated data is saved
SimulationDataFolder="/home/jose/MEGA/CPP_MCMC/Sparse-species-interactions/Macroecological-laws/Simulation-data"

##################################
#       FIXED PARAMETERS         #
##################################
#model

#Inverse growth rate, time to reach the steady state
tau <-  0.1
#Mean of the gaussian white noise
muNoise <-  0
#Standard deviation of the gaussian white noise
sigmaNoise <-  0.7
#Mean of the gaussian off diagonal term
muInt <-  0
#Standard deviation of the gaussian off diagonal term
sigmaInt <-  0.01
#Mean of the log-normal diagonal term
muLN <-  0.1
#Standard deviation of the log-normal diagonal term
sigmaLN <-  0.5

#Starting abundance of all species
abundanceStart <- 0.1

##################################
#       FIXED PARAMETERS         #
##################################

####################################################################
#             AFD: ABUNDANCE FLUCTUATION DISTRIBUTION              #
####################################################################

###########################################
#             S Experiment                #
###########################################

#Panel (A)

#Conectivty Fixed
C=0.5

#Sigma Noise Fixed
sigmaNoise=0.7

#Total number of species list
SList<-c(50,100,150)

#Read Data of the different experiments

#S experiment
AbundancesS=paste0("S_experiment/Abundances_C_",C,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)

#Complete Folders Name
AbundancesS_Folder=paste0(SimulationDataFolder,"/",AbundancesS)

#Generate and save data
expType="S"
AFD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (SList,expType,AbundancesS_Folder)

###########################################
#             S Experiment                #
###########################################

###########################################
#             C Experiment                #
###########################################

#Panel B

#Conectivty List
CList <- c(0.1,0.5,1.0)

#Sigma Noise Fixed
sigmaNoise=0.7

#Total number of species list
S <- 100

#Folder Name
AbundancesC=paste0("C_experiment/Abundances_S_",S,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)
AbundancesC_Folder=paste0(SimulationDataFolder,"/",AbundancesC)



#Generate and Save Data
expType="C"
AFD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (CList,expType,AbundancesC_Folder)

###########################################
#             C Experiment                #
###########################################

###########################################
#             Noise Experiment            #
###########################################

#Panel C

#Conectivty List
#C <- 0.5

#Sigma Noise Fixed
#NoiseList=c(0.1,0.25,0.5)

#Total number of species list
#S <- 100

#Folder Name
#AbundancesNoise=paste0("Noise_experiment/Abundances_S_",S,"_C_",C,"_tau_0.1_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)
#AbundancesNoise_Folder=paste0(SimulationDataFolder,"/",AbundancesNoise)


#Generate and Save data
#expType="sigmaNoise"
#AFD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (NoiseList,expType,AbundancesNoise_Folder)

###########################################
#             Noise Experiment            #
###########################################

###########################################
#             sigmaInt Experiment         #
###########################################

#Conectivty List
C <- 0.5

#Sigma Noise Fixed
sigmaNoise <- 0.7

#Total number of species list
S <- 100

#Interaction Strength List
sigmaInt_List <- c(0.01,0.02,0.03)

#Folder Name
Abundances_sigmaInt=paste0("sigmaInt_experiment/Abundances_S_",S,"_C_",C,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0")
Abundances_sigmaInt_Folder=paste0(SimulationDataFolder,"/",Abundances_sigmaInt)


#Generate and Save data
expType="sigmaInt"
AFD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (sigmaInt_List,expType,Abundances_sigmaInt_Folder)

###########################################
#             sigmaInt Experiment         #
###########################################

####################################################################
#             AFD: ABUNDANCE FLUCTUATION DISTRIBUTION              #
####################################################################

####################################################################
#             MAD: MEAN ABUNDANCE DISTRIBUTION                  #
####################################################################

###########################################
#             S Experiment                #
###########################################

#Panel (A)

#Conectivty Fixed
C=0.5

#Sigma Noise Fixed
sigmaNoise=0.7

#Total number of species list
SList<-c(50,100,150)

#Read Data of the different experiments

#S experiment
AbundancesS=paste0("S_experiment/Abundances_C_",C,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)

#Complete Folders Name
AbundancesS_Folder=paste0(SimulationDataFolder,"/",AbundancesS)



#Generate and save data
expType="S"
MAD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (SList,expType,AbundancesS_Folder)

###########################################
#             S Experiment                #
###########################################

###########################################
#             C Experiment                #
###########################################

#Panel B

#Conectivty List
CList <- c(0.1,0.5,1.0)

#Sigma Noise Fixed
sigmaNoise=0.7

#Total number of species list
S <- 100

#Folder Name
AbundancesC=paste0("C_experiment/Abundances_S_",S,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)
AbundancesC_Folder=paste0(SimulationDataFolder,"/",AbundancesC)



#Generate and Save Data
expType="C"
MAD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (CList,expType,AbundancesC_Folder)

###########################################
#             C Experiment                #
###########################################

###########################################
#             Noise Experiment            #
###########################################

#Panel C

#Conectivty List
#C <- 0.5

#Sigma Noise Fixed
#NoiseList=c(0.1,0.25,0.5)

#Total number of species list
#S <- 100

#Folder Name
#AbundancesNoise=paste0("Noise_experiment/Abundances_S_",S,"_C_",C,"_tau_0.1_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)
#AbundancesNoise_Folder=paste0(SimulationDataFolder,"/",AbundancesNoise)


#Generate and Save data
#expType="sigmaNoise"
#MAD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (NoiseList,expType,AbundancesNoise_Folder)

###########################################
#             Noise Experiment            #
###########################################


###########################################
#             sigmaInt Experiment         #
###########################################

#Conectivty List
C <- 0.5

#Sigma Noise Fixed
sigmaNoise <- 0.7

#Total number of species list
S <- 100

#Interaction Strength List
sigmaInt_List <- c(0.01,0.02,0.03)

#Folder Name
Abundances_sigmaInt=paste0("sigmaInt_experiment/Abundances_S_",S,"_C_",C,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0")
Abundances_sigmaInt_Folder=paste0(SimulationDataFolder,"/",Abundances_sigmaInt)


#Generate and Save data
expType="sigmaInt"
MAD_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (sigmaInt_List,expType,Abundances_sigmaInt_Folder)

###########################################
#             sigmaInt Experiment         #
###########################################

####################################################################
#             MAD: MEAN ABUNDANCE DISTRIBUTION                     #
####################################################################


####################################################################
#             TAYLOR: MEAN ABUNDANCE DISTRIBUTION                  #
####################################################################

###########################################
#             S Experiment                #
###########################################

#Panel (A)

#Conectivty Fixed
C=0.5

#Sigma Noise Fixed
sigmaNoise=0.7

#Total number of species list
SList<-c(50,100,150)

#Read Data of the different experiments

#S experiment
AbundancesS=paste0("S_experiment/Abundances_C_",C,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)

#Complete Folders Name
AbundancesS_Folder=paste0(SimulationDataFolder,"/",AbundancesS)



#Generate and save data
expType="S"
TAYLOR_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (SList,expType,AbundancesS_Folder)

###########################################
#             S Experiment                #
###########################################

###########################################
#             C Experiment                #
###########################################

#Panel B

#Conectivty List
CList <- c(0.1,0.5,1.0)

#Sigma Noise Fixed
sigmaNoise=0.7

#Total number of species list
S <- 100

#Folder Name
AbundancesC=paste0("C_experiment/Abundances_S_",S,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)
AbundancesC_Folder=paste0(SimulationDataFolder,"/",AbundancesC)



#Generate and Save Data
expType="C"
TAYLOR_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (CList,expType,AbundancesC_Folder)

###########################################
#             C Experiment                #
###########################################

###########################################
#             Noise Experiment            #
###########################################

#Panel C

#Conectivty List
#C <- 0.5

#Sigma Noise Fixed
#NoiseList=c(0.1,0.25,0.5)

#Total number of species list
#S <- 100

#Folder Name
#AbundancesNoise=paste0("Noise_experiment/Abundances_S_",S,"_C_",C,"_tau_0.1_muLN_0.1_sigmaLN_0.5_muInt_0_sigmaInt_",sigmaInt)
#AbundancesNoise_Folder=paste0(SimulationDataFolder,"/",AbundancesNoise)


#Generate and Save data
#expType="sigmaNoise"
#TAYLOR_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (NoiseList,expType,AbundancesNoise_Folder)

###########################################
#             Noise Experiment            #
###########################################

###########################################
#             sigmaInt Experiment         #
###########################################

#Conectivty List
C <- 0.1

#Sigma Noise Fixed
sigmaNoise <- 0.7

#Total number of species list
S <- 100

#Interaction Strength List
sigmaInt_List <- c(0.01,0.02,0.03)

#Folder Name
Abundances_sigmaInt=paste0("sigmaInt_experiment/Abundances_S_",S,"_C_",C,"_tau_0.1_sigmaNoise_",sigmaNoise,"_muLN_0.1_sigmaLN_0.5_muInt_0")
Abundances_sigmaInt_Folder=paste0(SimulationDataFolder,"/",Abundances_sigmaInt)


#Generate and Save data
expType="sigmaInt"
TAYLOR_SIMULATION_SINGLE_PLOT_ECOSYSTEMS_DATA (sigmaInt_List,expType,Abundances_sigmaInt_Folder)

###########################################
#             sigmaInt Experiment         #
###########################################

####################################################################
#             TAYLOR: MEAN ABUNDANCE DISTRIBUTION                  #
####################################################################