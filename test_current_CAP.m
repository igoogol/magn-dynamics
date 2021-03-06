clc;clear all
%% paramters
hbar = 6.582e-16;  % hbar, in unit of eV*s
e = 1.602e-19;     % charge of 1 electron, in unit of C
gamma_0 = 1.76e11; % gyromagnetic ratio of free electron, in unit of 1/(s*T)
Ms = 1.09e6;       % saturation magnitude, in unit of A/m
a = 0.045;         % lattice constant, in unit of nm 
t0 = 15.25;        % coupling between neighboring sites in central region
t1 = 1.0*t0;       % coupling between neighboring sites in leads
t2 = 0.8*t0;       % coupling between leads and central region
tSO = 0.098*t0;    % 0.4 eV - 3.3 eV (0.0262 t - 0.2164 t)
Ef = 0.484*t0;     % Fermi energy, 7.38 eV 
NB = 101;          % number of voltage points
NE = 4000;          % number of energy points
Bias = linspace(0.0,10.0*t0,NB); 
Energy = linspace(0.0001,9.0*t0,NE); 
kT = 0.001679*t0;  % temperature
M = 4;             % transverse layers
N = 5;             % translational layers 
K = 2;             % layers of thickness 
L = 20;             % no. of CAP layers in the leads
gamma_e = gamma_0*hbar/t0;    % gyromagnetic ratio (electrons), unit t0/(hbar*T)
gamma = gamma_e;              % gyromagnetic ratio (localized spin, or magnets)
V = (K-1)*(N-1)*(M-1)*(a*10^(-9))^3;  % volume of the system, unit m**3 
Sm = Ms*V/gamma_0/hbar/e;     % localized spin, dimensionless
J = 0.06557*t0;               % exchange interaction, ~ 1eV (0.066 t) 
JM = J*Sm;                    % magnetization 
D = 0.92*gamma_e/(2.0*Ms);    % anisotropic field parameter
TimeLength = 500000;          % Maximum total time scale
B = [-0.3,0.0,0.0]; % external field
theta = 5.0*pi/180.0;      % initial value of theta
phi = 0.0;                    % initial value of theta

NB = 1;
Bias = [0];

parameters_Hc = [t0,tSO,JM,M,N,gamma_e,B,D,Sm];
Hc = get_Hamiltonian_central(theta,phi,parameters_Hc);  

parameters_Hl = [t1,M,L];
Hl = get_Hamiltonian_lead(parameters_Hl);

%parameters_Ht = [t2,tSO,M];
parameters_Ht = [t2,0.0,M];
Ht = get_Hamiltonian_coupling(parameters_Ht);

%% Construct full Hamiltonian 
dimL = 2*M*L; % dimension of the leads (with CAP)
dimC = 2*M*N; % dimension of the central region
HLL = Hl;
HLC = [zeros(dimL-2*M,2*M),zeros(dimL-2*M,dimC-2*M);Ht,zeros(2*M,dimC-2*M)];
HLR = zeros(dimL,dimL);
HCL = HLC';
HCC = Hc;
HCR = [zeros(dimC-2*M,2*M),zeros(dimC-2*M,dimL-2*M);Ht',zeros(2*M,dimL-2*M)];
HRL = HLR';
HRC = HCR';
HRR = Hl;

H = [HLL,HLC,HLR;HCL,HCC,HCR;HRL,HRC,HRR];

% %% Get self-energy
% HL00 = [4.0*t1,0.0;0.0,4.0*t1];
% % HL00 = np.array([[4.0*t1,0.5*JM*np.exp(-1j*phi)],[0.5*JM*np.exp(1j*phi),4.0*t1]]) 
% HLT1 = [-t1,0.0;0.0,-t1];
% HLT2 = [-t1,0.0;0.0,-t1];
% HL0 = generate_block_tridiag(HL00,HLT1,M);
% HL1 = generate_block_diag(HLT2,M);
% Ht = [-t2,0.0;0.0,-t2];
% HT = generate_block_diag(Ht,M);
% 
% [SL,SR] = get_self_energy(Energy,M,HL0,HL1,HT);
      
%==========================================================================
%%                       Calculate current (CAP)
%==========================================================================
c = 2.62206;                   % CAP constant 
x = [1.12449e-1 0 8.28735e-3]; % CAP parameters
LL = (L-1)*a;                  % length of buffer layer
x1=a;   
x2=x1+LL;
% ys = linspace(0,0.98,L); 
% for ii = 1:L
%     y(2*M*(ii-1)+1:2*M*ii) = ys(ii);
% end
% fy=(4/(c^2)).*(-1./((1+y).^2)+1./((1-y).^2));
% for m=1:3
%     fy = fy + x(m)*(c^m)*(y.^m);
% end
% w=t0*(2*pi)^2*(a/LL)^2.*fy;

ys = linspace(0,0.99,L); 
for ii = 1:L
    y(2*M*(ii-1)+1:2*M*ii) = ys(ii);
end
fy=(4/(c^2)).*(1./((1+y).^2) + 1./((1-y).^2) -2); % y = z/deltaz
w=t0*(2*pi/(L-1))^2.*fy; % equivalent to w=t0*(2*pi)^2*(a/LL)^2.*fy;

% W matrix
wr = w;
wl = w(end:-1:1);
WR = diag(wr);
WL = diag(wl);

% Effective Hamiltonian
H = [HLL-1i*WL,HLC,HLR;HCL,HCC,HCR;HRL,HRC,HRR-1i*WR];

% WWL = zeros(length(H),length(H));
% WWL(1:length(WL),1:length(WL)) = WL;
% WWR = zeros(length(H),length(H));
% WWR(end-length(WL)+1:end,end-length(WL)+1:end) = WR;

dE = Energy(2) - Energy(1);
Current = [];
tic;
for jj = 1:NB
    Delta = Bias(jj);
    muL = 0.5*Delta; 
    muR = -0.5*Delta; 
    for ii = 1:NE
        E = Energy(ii);
        EI = E*eye(length(H));
        FL(ii) = cal_fermi(E,muL,Ef,kT);  % fermi function of left lead
        FR(ii) = cal_fermi(E,muR,Ef,kT);  % fermi function of right lead    
        Gr = inv(EI-H);                  % G retarded
        GLR = Gr(1:2*L*M,end-2*L*M+1:end);
        Trans(ii) = 4*trace(WL*GLR*WR*GLR');
        %Trans(ii) = 4*trace(WWL*Gr*WWR*Gr');
        clear EI Gr GLR
    end
%     Trans = funs.data_trans_trim(Trans,[])
    Ic = 0.0;
    for ii=1:NE
        Ic = Ic - (FL(ii)-FR(ii))*Trans(ii)*(dE/(2.0*pi));     % charge current 
    end
    fprintf('V = %f, Charge current is:%f \n',Delta,Ic);
    Current(jj) = Ic;
    
    %figure();
    plot(Energy,Trans);
    xlabel('Energy');
    ylabel('Transmission');
end
toc;
