
# Define parameters
H=100 # Physical length of the domain (100mts)
L=4 # Proportionality of the domain. number of columns(nx)=L*ny.
ny= 100 # Number of rows.
tol= 1.0e-6 # Poisson solver tolerance.
maxiter= 10 # Max num of iterations of the Poisson solver.
u_in= 7.0 # Inlet velocity. #
Re = 5000
rho = 1.22e-3 # Density of dry air at 1 atm and 15 C.
# rho = 1.0
tau = L * H / u_in # domain flow-through time
tf = 2.0 * tau # Simulation time