
# coding: utf-8

# ## Super Awesome Time evolution of 3-level System
# 
# Updated 12/3/2017

# In[24]:

# Import necessary modules & functions
# setup the matplotlib graphics library and configure it to show 
# figures inline in the notebook

get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from timeEvo_lambda import *


# In[25]:

#Python hdf5 library 

import h5py

f = h5py.File("timeEvo.hdf5", "a")


# In[140]:

# Parameters ('static' parameters are defined in timeEvo function)
# ================================================================

Omega = 1e-3*2*pi      # laser rabi frequency
 
gc = 1.e-3 * 2 * pi    # cavity coupling

n_th = 0.9    # avg number of thermal bath excitation
     
delta = 1.e-3 * 2 * pi  # detuning (exact value - not relative to Omega)

tlist = np.linspace(0,2e3,2000)   # total time to evolve system



# In[141]:

# Call timeEvo function and extract its output
# Usage: timeEvo(Omega, n_th, delta, gc , tlist, Hamiltonian, coupling)
# 0 for ladder system, 1 for lambda 
# 0 for cavity decoupled, 1 for coupled
# =================================================================== 

coupling = timeEvo(Omega, n_th, delta, gc, tlist, 1, 0)

output = timeEvo.output

n_c = output.expect[0]
n_g = output.expect[1]
n_r1 = output.expect[2]
n_r2 = output.expect[3]


# In[142]:

fig, axes = plt.subplots(1, 1, figsize=(10,6))
plt.gca().set_position((.15, .35, .8, .6)) # to make a bit of room for extra text

# Some variables for automated label generation
Delta = delta / (2 * pi)
tspan = max(tlist)
if coupling == 1:
    isCoupled = "On"
else:
    isCoupled = "Off"

# Graphs    
axes.plot(tlist, n_c, label="Cavity")
axes.plot(tlist, n_g, label="Population in |g>")
axes.plot(tlist, n_r1, label="Population in |r1>")
axes.plot(tlist, n_r2, label="Population in |r2>")
axes.legend(loc='upper right')
axes.set_xlabel('Time')
axes.set_ylabel('Occupation probability')
axes.set_title('Time evolution of system')
plt.figtext(.1, .15, "System = lambda \nOmega = 1 MHz \nDelta= %1.1e GHz \nCoupling = %1.0f \nNth= %1.2f" % (Delta, coupling,n_th) )
plt.ylim(ymin=-0.01)
plt.show()



# In[143]:

# Auto-creates groups and datasets in an HDF5 file.  Will error if new data would overwrite saved data.

dset_c = f.create_dataset('/lambda/nth=%1.2f/delta=%1.1eGHz/coupling%s/cav_population' % (n_th,Delta,isCoupled),data = n_c)
dset_g = f.create_dataset('/lambda/nth=%1.2f/delta=%1.1eGHz/coupling%s/g_population' % (n_th,Delta,isCoupled),data = n_g)
dset_r1 = f.create_dataset('/lambda/nth=%1.2f/delta=%1.1eGHz/coupling%s/r1_population' % (n_th,Delta,isCoupled),data = n_r1)
dset_r2 = f.create_dataset('/lambda/nth=%1.2f/delta=%1.1eGHz/coupling%s/r2_population' % (n_th,Delta,isCoupled),data = n_r2)

# A bunch of attributes which accompany the HDF5 data

f.attrs["omega"] = "1MHz"
f.attrs["wc"] = "15GHz"
f.attrs["wa"] = "15GHz"
f.attrs["Q"] = "1e5"
f.attrs["kappa"] = "150KHz"
f.attrs["gamma_r1"] = "3.4KHz"
f.attrs["gamma_r2"] = "1.5KHz"
f.attrs["fock_states"] = 60
f.attrs["initial_state"] = 3,0

dset_c.attrs["system"] = "lambdaSys"
dset_c.attrs["tspan"] = tspan
dset_c.attrs["n_th"] = "%1.2f" %n_th
dset_c.attrs["delta"] = "%1.1eGHz" %Delta
dset_c.attrs["coupling_str"] = "1MHz"
dset_c.attrs["coupling"] = "%s" %isCoupled

dset_g.attrs["system"] = "lambdaSys"
dset_g.attrs["tspan"] = tspan
dset_g.attrs["n_th"] = "%1.2f" %n_th
dset_g.attrs["delta"] = "%1.1eGHz" %Delta
dset_g.attrs["coupling_str"] = "1MHz"
dset_g.attrs["coupling"] = "%s" %isCoupled

dset_r1.attrs["system"] = "lambdaSys"
dset_r1.attrs["tspan"] = tspan
dset_r1.attrs["n_th"] = "%1.2f" %n_th
dset_r1.attrs["delta"] = "%1.1eGHz" %Delta
dset_r1.attrs["coupling_str"] = "1MHz"
dset_r1.attrs["coupling"] = "%s" %isCoupled

dset_r2.attrs["system"] = "lambdaSys"
dset_r2.attrs["tspan"] = tspan
dset_r2.attrs["n_th"] = "%1.2f" %n_th
dset_r2.attrs["delta"] = "%1.1eGHz" %Delta
dset_r2.attrs["coupling_str"] = "1MHz"
dset_r2.attrs["coupling"] = "%s" %isCoupled


# In[136]:

dset_c.attrs["delta"]


# In[137]:

dset_c.attrs["coupling"]


# In[138]:

dset_c.attrs["system"]


# In[139]:

dset_c.attrs["n_th"]

