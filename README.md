# Trabajo-Final-Senior-thesis
Repositorio con los códigos de Fortran que utilice para mi trabajo final y archivos .ipynb de Python que utilice para presentar los resultados obtenidos / Repository containing the Fortran codes I used for my final project and Python .ipynb files I used to present the obtained results.

random_module.f90: Module used to generate random numbers.

WS-bidimensional-euc-right.f90: Creates the Watts-Strogatz network in two dimensions with periodic boundary conditions and measures the free path and clustering coefficient of the network. This was done to check that the network possessed small-world properties.

GH_bid_tran.f90: Implements the Greenberg-Hastings dynamics on the Watts-Strogatz two-dimensional network and runs the dynamics for 1000 time steps, calculating the order parameter of the system to estimate relaxation times.

GH_bid.f90: Finally, implements the Greenberg-Hastings dynamics on the network and calculates the order parameter and susceptibility of the system.
