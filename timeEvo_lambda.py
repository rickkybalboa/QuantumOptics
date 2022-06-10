
# coding: utf-8

# ## Function for Time evolution of 3-level System, Ladder & Lambda
# 
# Updated 02/03/2017

import math
import time
from qutip import *
from math import *


def timeEvo(Omega, n_th, delta, gc, tlist, Hamiltonian,coupling):


    Omega = Omega
    n_th = n_th
    tlist = tlist
    
    delta = delta
    gc = gc
    coupling = coupling
    Hamiltonian = Hamiltonian
    
    if coupling == 1:
        gc = gc
    else:
        gc = 0

# Constant parameters

    wc = 15  * 2 * pi  # cavity frequency
    wa = 15  * 2 * pi  # atom frequency
    Q  = 1.e-4
    kappa = 1.5e-4       # cavity dissipation rate
    gamma = 0.1e-6 * 2 * pi       # atom dissipation rate
    gamma_r1 = 3.4e-6 * 2 *pi
    gamma_r2 = 1.5e-6 * 2 * pi
    N = 60              # number of cavity fock states
# intial state
    psi0 = tensor(basis(N,0), basis(3,0))    

    g  = tensor(qeye(N),basis(3,0))
    r1 = tensor(qeye(N),basis(3,1)) 
    r2 = tensor(qeye(N),basis(3,2))

# operators
    a  = tensor(destroy(N), qeye(3))
    sm = tensor(qeye(N), destroy(3))
    sig11 = g*g.dag()
    sig22 = r1*r1.dag()
    sig33 = r2*r2.dag()

# Ladder Hamiltonian
    H_c = delta * g*g.dag()+ Omega/2*(g*r1.dag()+r1*g.dag()) 
    H_jc = wc * a.dag() * a + wa * r2*r2.dag() + gc * (a * r2 * r1.dag() + a.dag() * r1 * r2.dag() )
    H_ladder = H_c + H_jc

# Lambda Hamiltonian
    H1 = delta*g*g.dag()+Omega/2*(r2*g.dag()+g*r2.dag())
    H2 = (wc * a.dag() * a) + (-wa * r1*r1.dag()) + gc * (a.dag() * (r1*r2.dag()) + a * (r2*r1.dag()))
    H_lambda = H1+H2
    
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

    rate = gamma_r1                                
    if rate > 0.0:
        c_ops.append(sqrt(rate) * g*r1.dag())   
    

    rate = gamma_r2                               
    if rate > 0.0:
        c_ops.append(sqrt(rate) * g*r2.dag())  
        
        t = time.time()

   # Master equation solver.  Chooses Hamiltonian to solve depending on input.
    rho0=tensor(thermal_dm(N,n_th),basis(3,0)*basis(3,0).dag())

    if Hamiltonian == 0:
        print("Solving Ladder System")
        timeEvo.output = mesolve(H_ladder, rho0, tlist, c_ops, [a.dag() * a, sig11, sig22, sig33])
    else:
        print("Solving Lambda System")
        timeEvo.output = mesolve(H_lambda, rho0, tlist, c_ops, [a.dag() * a, sig11, sig22, sig33])

    elapsed = time.time()-t
    print elapsed

    return coupling



