

# ## Steady State evolution of 3-level System
# 
# Updated 1/3/2017
# 
# Cavity coupling: 1MHz
# Omega=0.1MHz
# Coupling on

# setup the matplotlib graphics library and configure it to show 
# figures inline in the notebook
get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import numpy as np
import math
import time


# In[9]:

from qutip import *
from math import *


# In[10]:

wc = 15  * 2 * pi  # cavity frequency
wa = 15  * 2 * pi  # atom frequency
gc = 1.e-3 * 2 * pi  # coupling strength
Q  = 1.e-4
kappa = 1.5e-4       # cavity dissipation rate
gamma = 0.1e-6 * 2 * pi       # atom dissipation rate
gamma_r1 = 0.1e-6 * 2 *pi
gamma_r2 =0.1e-6 * 2 * pi
N = 60              # number of cavity fock states
Omega = 0.1e-3*2*pi 



# In[11]:



# intial state
psi0 = tensor(basis(N,0), basis(3,0))    # start with atom in ground state

g  = tensor(qeye(N),basis(3,0))
r1 = tensor(qeye(N),basis(3,1)) 
r2 = tensor(qeye(N),basis(3,2))

# operators
a  = tensor(destroy(N), qeye(3))
rho11 = g*g.dag()
rho22 = r1*r1.dag()
rho33 = r2*r2.dag()


# In[12]:

t = time.time()
#Define Calculation Range
detune=np.linspace(-20*Omega,20*Omega,250)
Nth=np.array([0,0.05,0.9,2.3,3.7,5.1])

# Initialise Output Arrays
A=np.zeros([Nth.size,detune.size])
B=np.zeros([Nth.size,detune.size])
n_c= np.zeros([Nth.size,detune.size])
n_r1=np.zeros([Nth.size,detune.size])
n_r2=np.zeros([Nth.size,detune.size])
n_g= np.zeros([Nth.size,detune.size])

for idx,n_th in enumerate(Nth):   

    print n_th
    

    # Assign Dissipation Terms  
    # ========================
    c_ops = []
    
    # cavity relaxation
    rate = kappa * (1 + n_th)
    if rate > 0.0:
        c_ops.append(math.sqrt(rate) * a)

    # cavity excitation, if temperature > 0
    rate = kappa * n_th
    if rate > 0.0:
        c_ops.append(math.sqrt(rate) * a.dag())

    # qubit relaxation
    #rate = gamma
    #if rate > 0.0:
        #    c_ops.append(math.sqrt(rate) * sm)
    
    rate = gamma_r1                                
    if rate > 0.0:
        c_ops.append(sqrt(rate) * g*r1.dag())   
    

    rate = gamma_r2                               
    if rate > 0.0:
        c_ops.append(sqrt(rate) * g*r2.dag())    
        
    #Determine Steadystate
    #rhoSS=steadystate(H,c_ops)
    #Extract steady-state cavity number
    #B[idx]=expect(rhoSS,a.dag()*a)    
        
        # Solve Steady State
    # ==================
    for ctr,o in enumerate(detune):
        #Define Hamiltonian
        H_c = o * g*g.dag()+ Omega/2*(g*r1.dag()+r1*g.dag()) 
        H_jc = wc * a.dag() * a + wa * r2*r2.dag() + gc * (a * r2 * r1.dag() + a.dag() * r1 * r2.dag() )
        H = H_c + H_jc        
        #Determine Steadystate
        rhoSS=steadystate(H,c_ops)
        #Extract steady-state cavity number
        A[idx,ctr]=expect(rhoSS,a.dag()*a)
        n_g[idx,ctr] = expect(rhoSS,g*g.dag()) 
        n_r1[idx,ctr] = expect(rhoSS,r1*r1.dag()) 
        n_r2[idx,ctr] = expect(rhoSS,r2*r2.dag())
    
elapsed = time.time()-t
print 'Time taken:',elapsed,'s'





