######################################################################
#                       LIBRARIES                                    #
######################################################################

import numpy as np
import os
from faux import InteractionMatrix_LNDiagonal, SLM_interactions_dynamics

#####################################################################
#                      LIBRARIES                                    #
#####################################################################

######################################################################
#               Total Numer of species (S) panel                     #
######################################################################

'''
Panel A) here I keep fixed the conectivity and  the noise intensity and
move the total number of species S
'''

#############################################
#		         PARAMETERS		             #
#############################################

#model

#Inverse growth rate, time to reach the steady state
tau = 0.1
#Mean of the gaussian white noise
muNoise = 0
#Standard deviation of the gaussian white noise
sigmaNoise = 0.7
#Mean of the gaussian off diagonal term
muInt = 0
#Standard deviation of the gaussian off diagonal term
sigmaInt = 0.01
#Mean of the log-normal diagonal term
muLN = 0.1
#Standard deviation of the log-normal diagonal term
sigmaLN = 1

#Starting abundance of all species
abundanceStart = 0.1

#List of total number of species of the synthetic ecosystem
SList = [50,100,150]
#List of the interaction network conectivity
C = 0.5

#total numer of simulations
nSim = 2000
#total number of interaction matrix 
nMat = 100

#time
tStart = 0
tEnd = 7
dt = 0.01
timeList = np.arange(tStart,tEnd+dt,dt)


#############################################
#		        PARAMETERS		            #
#############################################

#############################################
#               Output Folders              #
#############################################

#Main output folder
MainOutDirS = "/home/jose/MEGA/Figures/Simulation-data/S_experiment"

#If the folder don't exist, then create it
if os.path.isdir(MainOutDirS) == False:
    os.mkdir(MainOutDirS)

'''
Folders to save the abundance data
'''

#parameters output folder
ParamFoldS = "Abundances_"+"_".join(["C",str(C),"tau",str(tau),"sigmaNoise",str(sigmaNoise),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt),"sigmaInt",str(sigmaInt)])

MainParamOutDirS = MainOutDirS+"/"+ParamFoldS

#If the folder don't exist, then create it
if os.path.isdir(MainParamOutDirS) == False:
    os.mkdir(MainParamOutDirS)

'''
Folders to save the interaction Matrix
'''

ParamIntMatFoldS=  "IntMatrix_"+ "_".join(["C",str(C),"tau",str(tau),"sigmaNoise",str(sigmaNoise),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt),"sigmaInt",str(sigmaInt)])

MainParamIntMatOutDirS=MainOutDirS+"/"+ParamIntMatFoldS

#If the folder don't exist, then create it
if os.path.isdir(MainParamIntMatOutDirS) == False:
    os.mkdir(MainParamIntMatOutDirS)

#############################################
#               Output Folders              #
#############################################

#############################################
#               Simulation                  #
#############################################

#species loop
for S in SList:

    #abundance folder
    SFold = "/"+"_".join(["S",str(S)])
    SAbnFold = MainParamOutDirS+"/"+SFold
 	
    #If the folder don't exist, then create it
    if os.path.isdir(SAbnFold) == False:
        os.mkdir(SAbnFold)
    
    #Interaction Matrix folder
    SMatFold = MainParamIntMatOutDirS+"/"+SFold
    
    #If the folder don't exist, then create it
    if os.path.isdir(SMatFold) == False:
        os.mkdir(SMatFold)
    
    #matrix loop
    for i in range(0,nMat):
                
        print("Matrix %s with S=%s" % (i,S))
        
        #Abundance files
        AbnFile =  SAbnFold+"/"+"Abn_"+str(i)+".csv"
        
        #matrix files
        MatFile =  SMatFold+"/"+"MAT_"+str(i)+".csv"
        
        A = InteractionMatrix_LNDiagonal(S,C,muLN,sigmaLN,muInt,sigmaInt)
            
        #dynamics
        dynamics = SLM_interactions_dynamics (tau,sigmaNoise,muNoise,nSim,timeList,abundanceStart,A)
            
        crossCommunity = dynamics[-1,:,:]
            
        #save 1 
        np.savetxt(AbnFile,crossCommunity,delimiter=",")
        
        #save 2 
        np.savetxt(MatFile,A,delimiter=",")
    
    
#############################################
#               Simulation                  #
#############################################    

#######################################################################
#                Total Numer of species (S) panel                     #
#######################################################################

'''

'''

#######################################################################
##               Conectivity (C) panel                                #
#######################################################################     
      

'''
Panel B) here I keep fixed the total numer of species and  the noise intensity and
move the conectivity C
'''

#############################################
#		        PARAMETERS		            #
#############################################

#model

#Inverse growth rate, time to reach the steady state
tau = 0.1
#Mean of the gaussian white noise
muNoise = 0
#Standard deviation of the gaussian white noise
sigmaNoise = 0.7
#Mean of the gaussian off diagonal term
muInt = 0
#Standard deviation of the gaussian off diagonal term
sigmaInt = 0.01
#Mean of the log-normal diagonal term
muLN = 0.1
#Standard deviation of the log-normal diagonal term
sigmaLN = 0.5

#Starting abundance of all species
abundanceStart = 0.1

#Total number of species of the synthetic ecosystem
S = 100

#List of the interaction network conectivity
CList = [0.1,0.5,1]

#total numer of simulations
nSim = 2000

#total number of interaction matrix 
nMat = 50

#time
tStart = 0
tEnd = 7
dt = 0.01
timeList = np.arange(tStart,tEnd+dt,dt)


#############################################
#		        PARAMETERS		            #
#############################################

#############################################
#               Output Folders              #
#############################################

#Main output folder
MainOutDirC = "/home/jose/MEGA/Figures/Simulation-data/C_experiment"

#If the folder don't exist, then create it
if os.path.isdir(MainOutDirC) == False:
    os.mkdir(MainOutDirC)

'''
Folders to save the abundances
'''

#parameters output folder
ParamFoldC = "Abundances_"+ "_".join(["S",str(S),"tau",str(tau),"sigmaNoise",str(sigmaNoise),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt),"sigmaInt",str(sigmaInt)])

MainParamOutDirC = MainOutDirC+"/"+ParamFoldC

#If the folder don't exist, then create it
if os.path.isdir(MainParamOutDirC) == False:
    os.mkdir(MainParamOutDirC)


'''
Folders to save the interaction Matrix
'''

ParamIntMatFoldC=  "IntMatrix_"+ "_".join(["S",str(S),"tau",str(tau),"sigmaNoise",str(sigmaNoise),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt),"sigmaInt",str(sigmaInt)])

MainParamIntMatOutDirC=MainOutDirC+"/"+ParamIntMatFoldC

#If the folder don't exist, then create it
if os.path.isdir(MainParamIntMatOutDirC) == False:
    os.mkdir(MainParamIntMatOutDirC)

#############################################
#               Output Folders              #
#############################################

#############################################
#               Simulation                  #
#############################################

#Conectivity loop
for C in CList:

    #abundance folder
    CFold = "/"+"_".join(["C",str(C)])
    CAbnFold = MainParamOutDirC+"/"+CFold
 	
    #If the folder don't exist, then create it
    if os.path.isdir(CAbnFold) == False:
        os.mkdir(CAbnFold)
    
    #Interaction Matrix folder
    CMatFold = MainParamIntMatOutDirC+"/"+CFold
    
    #If the folder don't exist, then create it
    if os.path.isdir(CMatFold) == False:
        os.mkdir(CMatFold)
    
    #matrix loop
    for i in range(0,nMat):
                
        print("Matrix %s with C=%s" % (i,C))
        
        #Abundance files
        AbnFile =  CAbnFold+"/"+"Abn_"+str(i)+".csv"
        
        #matrix files
        MatFile =  CMatFold+"/"+"MAT_"+str(i)+".csv"
        
        A = InteractionMatrix_LNDiagonal(S,C,muLN,sigmaLN,muInt,sigmaInt)
            
        #dynamics
        dynamics = SLM_interactions_dynamics (tau,sigmaNoise,muNoise,nSim,timeList,abundanceStart,A)
            
        crossCommunity = dynamics[-1,:,:]
            
        #save 1 
        np.savetxt(AbnFile,crossCommunity,delimiter=",")
        
        #save 2 
        np.savetxt(MatFile,A,delimiter=",")
    
#############################################
#               Simulation                  #
#############################################    


#######################################################################
##               Conectivity (C) panel                                #
#######################################################################    

'''

'''

#######################################################################
##               Interactions (sigmaInt) panel                        #
#######################################################################     
       

'''
Panel C) here I keep fixed the total numer of species and  the conectivity and
move sigma of the interactions
'''

#############################################
#		        PARAMETERS		            #
#############################################

#model

#Inverse growth rate, time to reach the steady state
tau = 0.1
#Mean of the gaussian white noise
muNoise = 0
#Standard deviation of the gaussian white noise
sigmaNoise = 0.1
#Mean of the gaussian off diagonal term
muInt = 0
#Standard deviation of the gaussian off diagonal term
sigmaIntList = [0.01,0.02,0.03]
#Mean of the log-normal diagonal term
muLN = 0.1
#Standard deviation of the log-normal diagonal term
sigmaLN = 0.5

#Starting abundance of all species
abundanceStart = 0.1

#Total number of species of the synthetic ecosystem
S = 100

#List of the interaction network conectivity
C = 0.5

#total numer of simulations
nSim = 1000

#total number of interaction matrix 
nMat = 100

#time
tStart = 0
tEnd = 6
dt = 0.01
timeList = np.arange(tStart,tEnd+dt,dt)


#############################################
#		        PARAMETERS		            #
#############################################
        
#############################################
#               Output Folders              #
#############################################

#Main output folder
MainOutDirsigmaInt = "/home/jose/small-perturbation/Simulation-data/sigmaInt_experiment"

#If the folder don't exist, then create it
if os.path.isdir(MainOutDirsigmaInt) == False:
    os.mkdir(MainOutDirsigmaInt)

'''
Folders to save the Abundances
'''

#parameters output folder
ParamFoldsigmaInt = "Abundances_"+"_".join(["S",str(S),"C",str(C),"tau",str(tau),"sigmaNoise",str(sigmaNoise),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt)])

MainParamOutDirsigmaInt = MainOutDirsigmaInt+"/"+ParamFoldsigmaInt

#If the folder don't exist, then create it
if os.path.isdir(MainParamOutDirsigmaInt) == False:
    os.mkdir(MainParamOutDirsigmaInt)

'''
Folders to save the interaction Matrix
'''

ParamIntMatFoldsigmaInt=  "IntMatrix_"+ "_".join(["S",str(S),"C",str(C),"tau",str(tau),"sigmaNoise",str(sigmaNoise),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt)])

MainParamIntMatOutDirsigmaInt=MainOutDirsigmaInt+"/"+ParamIntMatFoldsigmaInt

#If the folder don't exist, then create it
if os.path.isdir(MainParamIntMatOutDirsigmaInt) == False:
    os.mkdir(MainParamIntMatOutDirsigmaInt)

#############################################
#               Output Folders              #
#############################################
            
#############################################
#               Simulation                  #
#############################################

#Interactions
for sigmaInt in sigmaIntList:

    #abundance folder
    sigmaIntFold = "_".join(["sigmaInt",str(sigmaInt)])
    sigmaIntAbnFold = MainParamOutDirsigmaInt+"/"+sigmaIntFold
 	
    #If the folder don't exist, then create it
    if os.path.isdir(sigmaIntAbnFold) == False:
        os.mkdir(sigmaIntAbnFold)
    
    #Interaction Matrix folder
    sigmaIntMatFold = MainParamIntMatOutDirsigmaInt+"/"+sigmaIntFold
    
    #If the folder don't exist, then create it
    if os.path.isdir(sigmaIntMatFold) == False:
        os.mkdir(sigmaIntMatFold)
    
    #matrix loop
    for i in range(0,nMat):
                
        print("Matrix %s with sigmaInt=%s" % (i,sigmaInt))
        
        #Abundance files
        AbnFile =  sigmaIntAbnFold+"/"+"Abn_"+str(i)+".csv"
        
        #matrix files
        MatFile =  sigmaIntMatFold+"/"+"MAT_"+str(i)+".csv"
        
        A = InteractionMatrix_LNDiagonal(S,C,muLN,sigmaLN,muInt,sigmaInt)
            
        #dynamics
        dynamics = SLM_interactions_dynamics (tau,sigmaNoise,muNoise,nSim,timeList,abundanceStart,A)
            
        crossCommunity = dynamics[-1,:,:]
            
        #save 1 
        np.savetxt(AbnFile,crossCommunity,delimiter=",")
        
        #save 2 
        np.savetxt(MatFile,A,delimiter=",")
    
#############################################
#               Simulation                  #
#############################################   
        
######################################################################
#               Interactions (sigmaInt) panel                        #
######################################################################     
      

#######################################################################
##               Sigma Noise (sigmaNoise) panel                       #
#######################################################################     
      

'''
Panel C) here I keep fixed the total numer of species and  the conectivity and
move sigma intensity noise
'''

#############################################
#		        PARAMETERS		            #
#############################################

#model

#Inverse growth rate, time to reach the steady state
tau = 0.1
#Mean of the gaussian white noise
muNoise = 0
#Standard deviation of the gaussian white noise
sigmaNoiseList = [0.2,0.4,0.8]
#Mean of the gaussian off diagonal term
muInt = 0
#Standard deviation of the gaussian off diagonal term
sigmaInt = 0.03
#Mean of the log-normal diagonal term
muLN = 0.1
#Standard deviation of the log-normal diagonal term
sigmaLN = 0.5

#Starting abundance of all species
abundanceStart = 0.1

#Total number of species of the synthetic ecosystem
S = 100

#List of the interaction network conectivity
C = 0.5

#total numer of simulations
nSim = 2000

#total number of interaction matrix 
nMat = 50

#time
tStart = 0
tEnd = 7
dt = 0.01
timeList = np.arange(tStart,tEnd+dt,dt)


#############################################
#		        PARAMETERS		            #
#############################################
        
#############################################
#               Output Folders              #
#############################################

#Main output folder
MainOutDirNoise = "/home/jose/MEGA/Figures/Simulation-data/Noise_experiment"

#If the folder don't exist, then create it
if os.path.isdir(MainOutDirNoise) == False:
    os.mkdir(MainOutDirNoise)

'''
Folders to save the Abundances
'''

#parameters output folder
ParamFoldNoise = "Abundances_"+"_".join(["S",str(S),"C",str(C),"tau",str(tau),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt),"sigmaInt",str(sigmaInt)])

MainParamOutDirNoise = MainOutDirNoise+"/"+ParamFoldNoise

#If the folder don't exist, then create it
if os.path.isdir(MainParamOutDirNoise) == False:
    os.mkdir(MainParamOutDirNoise)

'''
Folders to save the interaction Matrix
'''

ParamIntMatFoldNoise=  "IntMatrix_"+ "_".join(["S",str(S),"C",str(C),"tau",str(tau),"muLN",str(muLN),"sigmaLN",str(sigmaLN),"muInt",str(muInt),"sigmaInt",str(sigmaInt)])

MainParamIntMatOutDirNoise=MainOutDirNoise+"/"+ParamIntMatFoldNoise

#If the folder don't exist, then create it
if os.path.isdir(MainParamIntMatOutDirNoise) == False:
    os.mkdir(MainParamIntMatOutDirNoise)

#############################################
#               Output Folders              #
#############################################
            
#############################################
#               Simulation                  #
#############################################

#Noise loop
for sigmaNoise in sigmaNoiseList:

    #abundance folder
    NoiseFold = "_".join(["sigmaNoise",str(sigmaNoise)])
    NoiseAbnFold = MainParamOutDirNoise+"/"+NoiseFold
 	
    #If the folder don't exist, then create it
    if os.path.isdir(NoiseAbnFold) == False:
        os.mkdir(NoiseAbnFold)
    
    #Interaction Matrix folder
    NoiseMatFold = MainParamIntMatOutDirNoise+"/"+NoiseFold
    
    #If the folder don't exist, then create it
    if os.path.isdir(NoiseMatFold) == False:
        os.mkdir(NoiseMatFold)
    
    #matrix loop
    for i in range(0,nMat):
                
        print("Matrix %s with sigmaNoise=%s" % (i,sigmaNoise))
        
        #Abundance files
        AbnFile =  NoiseAbnFold+"/"+"Abn_"+str(i)+".csv"
        
        #matrix files
        MatFile =  NoiseMatFold+"/"+"MAT_"+str(i)+".csv"
        
        A = InteractionMatrix_LNDiagonal(S,C,muLN,sigmaLN,muInt,sigmaInt)
            
        #dynamics
        dynamics = SLM_interactions_dynamics (tau,sigmaNoise,muNoise,nSim,timeList,abundanceStart,A)
            
        crossCommunity = dynamics[-1,:,:]
            
        #save 1 
        np.savetxt(AbnFile,crossCommunity,delimiter=",")
        
        #save 2 
        np.savetxt(MatFile,A,delimiter=",")
    
#############################################
#               Simulation                  #
#############################################   
        
#######################################################################
##               Sigma Noise (sigmaNoise) panel                       #
#######################################################################           
        
        





