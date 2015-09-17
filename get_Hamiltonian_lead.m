function [Hl] = get_Hamiltonian_lead(parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get Hamiltonian of the isolated leads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = parameters(1);
M = parameters(2);   % no. of tranverse modes in the leads
L = parameters(3);   % no. of (CAP) layers in the leads
HL00 = [4.0*t1,0.0;0.0,4.0*t1];
HLT1 = [-t1,0.0;0.0,-t1];
HLT2 = [-t1,0.0;0.0,-t1];
HL0 = generate_block_tridiag(HL00,HLT1,M);
HL1 = generate_block_diag(HLT2,M);
Hl = generate_block_tridiag(HL0,HL1,L);
