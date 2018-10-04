# Deriving_mesoscopic_equations

This repositiory contains codes used for simulating the pairwise interaction and the ternary
interaction model of collective behavior that are part of the book chapter entitiled
"Deriving mesoscopic models of collective behaviour for finite populations", to be
published in the Elsevier: Handbook of Statistics, Vol 40.

All codes run on MATLAB.

SSA_pairwise.m simulates the pairwise interaction model of collective behaviour and SSA_ternary.m simulates the ternary interaction model of collective behaviour using the the stochastic simulation algorithm. Both these codes run simulations for the time steps given by the user. The collective state is quantified as group polarization at constant epochs of time. The resulting data is a time series of group polarization. The codes also print the distribution of polarization.

SDE_simulate_pairwise.m and SED_simulate_ternary.m perform numerical integration of the mesoscopic stochastic differential equation describing the dynamics of group polarization for the pairwise and the ternary interaction model of collective behaviour. The output is the time series data of group polarization and the distribution of the same.

