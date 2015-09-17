function [Ht] = get_Hamiltonian_coupling(parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get Hamiltonian of the isolated leads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2 = parameters(1);
tSO = parameters(2);
M = parameters(3); % no. of transverse sites
H0 = [-t2,tSO;-tSO,-t2];
Ht = generate_block_diag(H0,M);