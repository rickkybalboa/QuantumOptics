
# # Fink Paper Modelling
# 
# An attempt to replicate the computational results obtained in 'Quantum-To-Classical Transition in Cavity Quantum Electrodynamics', Fink et al, 2010

# ## Transmission values with parameters in Fink paper
# 
# Here the transmission is found for several values of $n_{th}$ using the parameters given in the Fink paper.  The transmission is scaled over the max transmission value for an empty cavity


# Cavity Parameters (From FINK)
# =======================
N = 70                    # number of cavity fock states

wc = 6.44  * 2 * pi       # cavity frequency
wa = 6.44  * 2 * pi       # atom frequency
g  = 54.e-3 * 2 * pi      # coupling strength  ~54MHz
kappa = 3.2e-3 * 2 * pi   # cavity dissipation rate ~3.2MHz
gamma = 0.6e-3 * 2 * pi   # atom dissipation rate ~0.6MHz
eta = kappa*sqrt(0.01) * 2 * pi       # probe tone frequency ~1MHz
omega=np.linspace(6.35,6.55,250)* 2 * pi
Nth=np.array([0.05,0.09,0.19,0.42,1.1,2.7,8.4])


t=time.time()

# Define Operators
# ================
a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))

#Define Calculation Range
omega=np.linspace(6.35,6.55,101)* 2 * pi
Nth=np.array([0.05,0.09,0.19,0.42,1.1,2.7,8.4])

# Initialise Output Array
A=np.zeros([Nth.size,omega.size])
# Initialise Output Array
B=np.zeros([Nth.size])

for idx,n_th in enumerate(Nth):   

    print n_th
    # Assign Dissipation Terms  
    # ========================
    c_ops = []
    
    #cavity relaxation
    rate = kappa * (1 + n_th)
    if rate > 0.0:
        c_ops.append(sqrt(rate) * a)
        
    # cavity excitation, if temperature > 0
    rate = kappa * n_th
    if rate > 0.0:
        c_ops.append(sqrt(rate) * a.dag())
        
    #qubit relaxation    
    rate = gamma
    if rate > 0.0:
        c_ops.append(sqrt(rate) * sm) 
        
    #Define Hamiltonian
    H = eta*(a+a.dag())
    #Determine Steadystate
    rhoSS=steadystate(H,c_ops)
    #Extract steady-state cavity number
    B[idx]=expect(rhoSS,a.dag()*a)    

    # Solve Steady State
    # ==================
    for ctr,o in enumerate(omega):
        #Define Hamiltonian
        H = (wc-o) * a.dag() * a + (wa-o) * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())+eta*(a+a.dag())
        #Determine Steadystate
        rhoSS=steadystate(H,c_ops)
        #Extract steady-state cavity number
        A[idx,ctr]=expect(rhoSS,a.dag()*a)
    
elapsed = time.time()-t
print 'Time taken:',elapsed
