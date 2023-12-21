# In[1]:

#LIBRARIES


import numpy as np
import math as mt
import scipy.optimize as optimization
from scipy import special as sp
from scipy import stats

# In[2]:


def InteractionMatrix_LNDiagonal (S,C,muLN,sigmaLN,muInt,sigmaInt):
    
    """
    Interaction matrix with LogNormal (LN) diagonal and Gaussian off-diagonal fluctuations;
    
    INPUT
    - S (float): species number;
    - C (float): connectance;
    - muLN (float): LN mean of the diagonal elements; 
    - sigmaLN (float): LN widths of the diagonal elements;
    - muInt (float): Gaussian mean of the off-diagonal elements;
    - sigmaInt (float): Gaussian variance of the off-diagonal elements;
    
    OUTPUT
    - M (SxS array): i,j element provides the interaction weight between the species i and j
    """
    
    #interaction
    A = np.random.normal(muInt,sigmaInt,size=(S,S))
    
    #connectance
    boolMat = np.random.uniform(0,1,size=(S,S))
    A = (boolMat<C).astype(int)*A
    A = A-np.diag(np.diag(A))
    
    #diagonal
    K = 1/np.random.lognormal(muLN,sigmaLN,size=S)
    M = A-np.diag(K)
    
    return(M)


def LogisticFun (t,tau,K,x0):
    
    """
    Analytical solution of the logistic equation (without noise);
    
    INPUT
    - t (float): time instant;
    - tau (float): intrinsic growth time;
    - K (float): carrying capacity;
    - x0 (float): abundance at t=0;
    
    OUTPUT
    - x (float): abundance at time t;
    """
    
    x = (K*x0*mt.exp(t/tau))/((K-x0)+x0*mt.exp(t/tau))
    
    return x


# In[4]:


def SLM_dynamics (K,tau,sigmaNoise,muNoise,nSim,timeList,abundanceStart,S):
    
    """
    Stochastic Logistic equation (SLM) solution by Euler-Maruyama method;
    
    INPUT
    - K (float): carrying capacitities;
    - tau (float): intrinsic growth time;
    - sigmaNoise (float): variance of the Gaussian noise;
    - muNoise (float): mean value of the Gaussian noise;
    - nSim (float): number of simulation;
    - timeList (1d array): list of T time instants;
    - abundanceStart (float): initial abundance;
    - S (float): number of species;
    
    OUTPUT
    - longitudinal (T x S array): abundance of the species j at time i;
    - crossSectional (S x nSim array): abundance of the k trajectory at the last time for the j trajectory;
    """
  

    #initial state
    x0 = np.full(S,abundanceStart)

    #simulation loop

    abundanceMat = np.zeros((timeList.size,S,nSim),dtype=float)
    
    dt = round((timeList[-1]-timeList[0])/(len(timeList)),2)
    
    for sim in range(nSim):
        
        abundanceMat[0,:,sim] = x0
        
        #time loop
        for t in range(1,timeList.size):

            #dynamics
            x = abundanceMat[t-1,:,sim]
            noise = np.random.normal(muNoise,sigmaNoise,size=(S)) #Gaussian noise
            xt = x+x*dt/tau-x**2/(tau*K)*dt+(noise*x*mt.sqrt(dt)) #equation
        
            abundanceMat[t,:,sim]=np.round(xt,5)
    
    longitudinal = abundanceMat[:,:,0]
    crossSectional = abundanceMat[-1,:,:]
    
    return([longitudinal,crossSectional])



def SLM_dynamics_sigmaDist (K,tau,sigma,nSim,sigmaNoise,muNoise,timeList,abundanceStart,S):
    
    """
    Stochastic Logistic equation (SLM) solution by Euler-Maruyama method;
    
    INPUT
    - K (float): carrying capacitities;
    - tau (float): intrinsic growth time;
    - sigma (float): array of noise parameters
    - nSim (float): number of simulation;
    - timeList (1d array): list of T time instants;
    - abundanceStart (float): initial abundance;
    - S (float): number of species;
    
    
    OUTPUT
    - longitudinal (T x S array): abundance of the species j at time i;
    - crossSectional (S x nSim array): abundance of the k trajectory at the last time for the j trajectory;
    """
  
    

    #initial state
    x0 = np.full(S,abundanceStart)

    #simulation loop

    abundanceMat = np.zeros((timeList.size,S,nSim),dtype=float)
    
    dt = round((timeList[-1]-timeList[0])/(len(timeList)),2)
    
    for sim in range(nSim):
        
        abundanceMat[0,:,sim] = x0
        
        #time loop
        for t in range(1,timeList.size):

            #dynamics
            x = abundanceMat[t-1,:,sim]
            noise = np.random.normal(muNoise,sigmaNoise,size=(S)) #Gaussian noise
            xt = x+x*dt/tau-x**2/(tau*K)*dt+(noise*x*mt.sqrt(dt)) #equation
        
            abundanceMat[t,:,sim]=np.round(xt,5)
    
    longitudinal = abundanceMat[:,:,0]
    crossSectional = abundanceMat[-1,:,:]
    
    return([longitudinal,crossSectional])

# In[5]:


def LV_SLM_dynamics (K,tau,sigmaNoise,muNoise,nSim,timeList,abundanceStart,muInt,sigmaInt,S,C):
    
    """
    Stochastic Logistic equation (SLM) with Lotka-Volterra interaction (LV) solved by Euler-Maruyama method;
    
    INPUT
    - K (float): carrying capacitities;
    - tau (float): intrinsic growth time;
    - sigmaNoise (float): variance of the Gaussian noise;
    - muNoise (float): mean value of the Gaussian noise;
    - nSim (float): number of simulation;
    - timeList (1d array): list of T time instants;
    - abundanceStart (float): initial abundance;
    - muInt (float): mean value of the interaction weights;
    - sigmaInt (float): variance of the interaction weights;
    - S (float): number of species;
    - C (float): interaction network connectance;
    
    OUTPUT
    - longitudinal (T x S array): abundance of the species j at time i;
    - crossSectional (T x S x nSim array): abundance of the k trajectory at time i for the j trajectory;
    - steadyState (1d array of S-dim): stationary abundance calculated analytically;
    """

    #initial state
    x0 = np.full(S,abundanceStart)

    #network
    A = np.random.normal(muInt,sigmaInt,size=(S,S))
    boolMat = np.random.uniform(0,1,size=(S,S))
    A = (boolMat<C).astype(int)*A
    A = A-np.diag(np.diag(A))
    
    #stationary state
    steadyState = K*np.matmul(np.linalg.inv(np.identity(S)-K*A),np.ones(S))

    #simulation

    abundanceMat = np.zeros((timeList.size,S,nSim),dtype=float)
    
    dt = round((timeList[-1]-timeList[0])/(len(timeList)),2)
    
    #simulation loop
    for sim in range(nSim):
        
        abundanceMat[0,:,sim] = x0
        
        #time loop
        for t in range(1,timeList.size):
            
            #dynamics
            x = abundanceMat[t-1,:,sim]
            noise = np.random.normal(muNoise,sigmaNoise,size=(S)) #noise
            xt = x+x*dt/tau+(noise*x*mt.sqrt(dt))-x**2/(tau*K)*dt+x*np.matmul(A,x)*dt/tau #equation
            
            abundanceMat[t,:,sim]=xt
    
    longitudinal = abundanceMat[:,:,0]
    crossSectional = abundanceMat
    
    return([longitudinal,crossSectional,steadyState])
    


# In[6]:


def communityLongitudinal (timeList, abundanceMat, tSt, nComm):
    
    """
    Community-Species matrix as introduced by J. Grilli (2020)
    
    INPUT
    - timeList (1d array): list of T time instants;
    - abundanceMat (T x S array): mean abundance of the species j at time i averaged over all the simulations;
    - tSt (float): time to match the steady state;
    - nComm (float): community number
    
    
    OUTPUT
    - commSpeciesMat (S x nComm): abundance of i species in the j community
    """
    
    
    #stationary regime
    timeStList = timeList[timeList>tSt]
    
    #first stationary instant
    timeStNum = len(timeStList)
    
    #community index list
    commIndexList = np.linspace(len(timeList)-timeStNum,len(timeList)-1,nComm,dtype=int)
    
    #community-species matrix
    commSpeciesMat = np.transpose(np.take(abundanceMat,commIndexList,axis=0))
    
    return(commSpeciesMat)


# In[ ]:


def SLM_interactions_dynamics (tau,sigmaNoise,muNoise,nSim,timeList,abundanceStart,A):
    
    """
    Stochastic Logistic equation (SLM) with generic interaction solved by Euler-Maruyama method;
    
    INPUT
    - tau (float): intrinsic growth time;
    - sigmaNoise (float): variance of the Gaussian noise;
    - muNoise (float): mean value of the Gaussian noise;
    - nSim (float): number of simulation;
    - timeList (1d array): list of T time instants;
    - abundanceStart (float): initial abundance;
    - A (SxS array): generic interaction matrix
    
    OUTPUT
    - abundanceMat (T x S x nSim array): abundance of the species j at time i of k trajectory;
    """
    
    #DYNAMICS
    
    #species number
    S = A.shape[0]
    
    #initial state
    x0 = np.full(S,abundanceStart)

    #simulation

    abundanceMat = np.zeros((timeList.size,S,nSim),dtype=float)
    
    dt = round((timeList[-1]-timeList[0])/(len(timeList)),2)
    
    #simulation loop
    for sim in range(nSim):
        
        abundanceMat[0,:,sim] = x0
        
        #time loop
        for t in range(1,timeList.size):
            
            #dynamics
            x = abundanceMat[t-1,:,sim]
            noise = np.random.normal(muNoise,sigmaNoise,size=(S)) #noise
            xt = x+x*dt/tau+(noise*x*mt.sqrt(dt))+x*np.matmul(A,x)*dt/tau #equation
            
            abundanceMat[t,:,sim]=xt
    
    return(abundanceMat)
    


###################
def gamma_fit(x, b, k):
            # log10( sqrt(trigamma(b))/gamma(b)*(k)^b*exp( (b)*( x*sqrt( trigamma(b) )+digamma(b)+log(k) )-k*exp( ( x*sqrt( trigamma(b) )+digamma(b)+log(k) ) ) ) ) }
    return np.log10( np.sqrt(sp.polygamma(1,b))/sp.gamma(b)*((k)**b)*np.exp( (b)*( x*np.sqrt( sp.polygamma(1,b) )+sp.polygamma(0,b)+np.log(k) )-k*np.exp( ( x*np.sqrt( sp.polygamma(1,b) )+sp.polygamma(0,b)+np.log(k) ) ) ) ) 
    
def linear_fit(x,m,n):
    
    return m*x+n

def AFD_from_AbundanceTable(crossCommunity,bins_list):
    
    '''
    Function to compute the AFD from a abundance table
    -crossCommunity: abundance table, rows are species, columns are samples
    -bins_List: list of bins
    '''
    lp=np.log(crossCommunity)
    species_log_mean=np.mean(lp,axis=1)
    species_log_std=np.std(lp,axis=1)
    z= (lp-species_log_mean[:,None]) / species_log_std[:,None] 

    z=z.flatten()
    z=z[z>-4]

    AFD,bin_edges=np.histogram(z,bins_list,density=True)
      
    #popt, pcov=optimization.curve_fit(gamma_fit, bin_edges[0:nbins], np.log10(AFD), (1,1))
   
    return AFD
    

def MAD_from_AbundanceTable(crossCommunity,bins_list):
    
    '''
    Function to compute the MAD from a abundance table
    -crossCommunity: abundance table, rows are species, columns are samples
    -bins_List: list of bins
    '''

    lp=np.log(crossCommunity)
    species_log_mean=np.mean(lp,axis=1)
    z=(species_log_mean-np.mean(species_log_mean))/np.std(species_log_mean)
    z=z.flatten()
    z=z[z>-4]
    MAD,bin_edges=np.histogram(z,bins_list,density=True)
    
    return MAD
	

###############################################################################

def TAYLOR_from_AbundanceTable(crossCommunity):
    
    '''
    Function to compute the MAD from a abundance table
    -crossCommunity: abundance table, rows are species, columns are samples
    '''
    
    species_mean=np.mean(crossCommunity,axis=1)
    species_var=np.var(crossCommunity,axis=1)
    
    return species_mean[(species_mean>1e-6) & (species_var>1e-6)],species_var[(species_mean>1e-6) & (species_var>1e-6)]
    #return species_mean,species_var

# In[7]:


#def SLM_interactions_dynamics (tau,sigmaNoise,muNoise,nSim,nComm,timeList,abundanceStart,A):
    
    #"""
    ##Stochastic Logistic equation (SLM) with generic interaction solved by Euler-Maruyama method;
    
    ##INPUT
    #- tau (float): intrinsic growth time;
    #- sigmaNoise (float): variance of the Gaussian noise;
    #- muNoise (float): mean value of the Gaussian noise;
    #- nSim (float): number of simulation;
    #- nComm (float): number of simulation;
    #- timeList (1d array): list of T time instants;
    #- abundanceStart (float): initial abundance;
    #- A (SxS array): generic interaction matrix
    
    ##OUTPUT
    #- longitudinal (T x S array): abundance of the species j at time i (equispaced randomly selected);
    #- crossSectional (S x nSim array): abundance of the k trajectory at the last time for the j trajectory;
    #- steadyState (1d array of S-dim): stationary abundance calculated analytically;
    #"""
    
    ##DYNAMICS
    
    ##species number
    #S = A.shape[0]
    
    ##initial state
    #x0 = np.full(S,abundanceStart)

    ##stationary state
    #steadyState = -np.matmul(np.linalg.inv(A),np.ones(S))

    ##simulation

    #abundanceMat = np.zeros((timeList.size,S,nSim),dtype=float)
    
    #dt = round((timeList[-1]-timeList[0])/(len(timeList)),2)
    
    #simulation loop
    #for sim in range(nSim):
        
        #abundanceMat[0,:,sim] = x0
        
        #time loop
        #for t in range(1,timeList.size):
            
            #dynamics
            #x = abundanceMat[t-1,:,sim]
            #noise = np.random.normal(muNoise,sigmaNoise,size=(S)) #noise
            #xt = x+x*dt/tau+(noise*x*mt.sqrt(dt))+x*np.matmul(A,x)*dt/tau #equation
            
            #abundanceMat[t,:,sim]=xt
    
    
    ##COMMUNITY MATRIIX
    
    ##crossectional
    #crossSectional = abundanceMat[-1,:,:]
    
    ##longitudinal
    
    ##sort species according stationary abundances
    #arr1inds = steadyState.argsort()
    #steadyState = steadyState[arr1inds[::-1]]
    #abundanceMat = abundanceMat[:,arr1inds[::-1]]
    
    ##average over all simulations
    #abundanceAv = abundanceMat.mean(axis=2)
    
    ##equilibrium time
    #xAn = steadyState[0] 
    #tSt = timeList[abs(abundanceAv[:,0]-xAn)/xAn<1/100][0]
    #select the most abundant species: it is the last to reach equilibrium
    
    ##longitudina community matrix
    #longitudinal = communityLongitudinal (timeList, abundanceMat[:,:,0], tSt, nComm)
    
    #return([longitudinal,crossSectional,steadyState])
    


# In[ ]:




