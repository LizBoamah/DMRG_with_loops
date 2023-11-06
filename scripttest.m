clear
clc
% Example script to run the two_site_dmrg function

% Parameters for the Hubbard model
N = 6;             % Number of sites
P = 2;             % Local space dimension (spin up and spin down)
D = 10;            % Bond dimension
U = 0;             % Onsite interaction strength
t = 1;             % Hopping parameter
max_sweeps = 20;   % Maximum number of sweeps
mu = 0;
 dir = 'mixed';
% Run the two_site_dmrg function
[dmrg,E_exact, energy_values,new_energy_values, E_DMRG_first,E_DMRG_second] = two_site_dmrg(N, mu,P, D,t, max_sweeps);



