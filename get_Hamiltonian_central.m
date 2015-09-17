function [Hc] = get_Hamiltonian_central(theta,phi,parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get Hamiltonian of the central region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = parameters(1);
tSO = parameters(2);
JM = parameters(3);
M = parameters(4);
N = parameters(5);  
gamma_e = parameters(6);
Bx = parameters(7);
By = parameters(8);
Bz = parameters(9);
D = parameters(10);
Sm = parameters(11);
Smz = Sm*cos(theta);
H00 = [4.0*t0+0.5*JM*cos(theta), ...
       0.5*JM*sin(theta)*exp(-1j*phi); ...
       0.5*JM*sin(theta)*exp(1j*phi), ...
       4.0*t0-0.5*JM*cos(theta)] ...
     +gamma_e*[Bz, Bx-1i*By; Bx+1i*By, -Bz] ...
     +[-D*Smz^2, 0.0; 0.0, -D*Smz^2];
HT1 = [-t0,-1i*tSO;-1i*tSO,-t0];
HT2 = [-t0,tSO;-tSO,-t0];
Hc = get_Hamiltonian(H00,HT1,HT2,M,N);










