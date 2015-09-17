function [is_div,grR] = cal_surface_green_R(E,H0,H1) 
%==========================================================================
%  Calculate the surface Green's function of right lead: gr00. 
%  The algorithm is taken from Qiao.Z.H's PhD thesis P41-42. 
%  where:
%  err: accuracy
%  tt0: tilde{t}_0
%  TT0 = tt0*tt1*tt2*...*tt(n-1)  
%==========================================================================

% calculate right:
is_div = 1;
err = 1e-8;
dim = length(H0);
I = eye(dim);
eta = 1e-8;
while is_div == 1
    EI = (E+eta*1i)*eye(dim);
    t0 = inv(EI-H0)*H1';
    tt0 = inv(EI-H0)*H1;
    TR = tt0; % Don't assign TR = tt0 directly
    T = t0;   % Don't assign T = t0 directly
    t = eye(dim);     
    while sum(abs(t)) > err
        T0 = T;
        s1 = I - t0*tt0 - tt0*t0;
        s1_inv = inv(s1);
        s2 = t0*t0;
        s3 = tt0*tt0;
        t = s1_inv*s2;
        tt = s1_inv*s3;
        T = T0 + TR*t;  % update T
        t0 = t; % Don't assign t0 = t directly
        tt0 = tt; % Don't assign tt0 = tt directly
        TR = TR*tt;
        if isnan(sum(T))
            fprintf('T reaches singularity. Energy = %f, eta = %f',E,eta);
            eta = eta*1.05;  % update value of eta
            break
        elseif eta>1e-5
            fprintf('eta exceeds 1e-5. Energy = %f, eta = %f',E,eta);
            T = eye(dim);
            is_div = 0; % calculate grR incorrectly
        else
            is_div = 0; % calculate grR correctly 
        end
    end
end
grR = inv(EI-H0-H1*T);
clear err dim I EI t0 tt0 TR T T0 s1 s1_inv s2 s3