
# Define parameters
H=1 # Physical length of the domain (100mts)
L=4 # Proportionality of the domain. number of columns(nx)=L*ny.
tf = 10 # Simulation time.  
ny= 100 # Number of rows. 
tol= 1.0e-6 # Poisson solver tolerance. 
maxiter= 100 # Max num of iterations of the Poisson solver. 
u_in= 1.0 # Inlet velocity. 
nu=0.1
F=-60