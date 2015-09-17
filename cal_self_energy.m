function [is_div_L,is_div_R,SigmaL,SigmaR] = cal_self_energy(E,H0,H1,HT)
%==========================================================================
%  Calculate the surface Green's function of right lead: gr00. 
%  The algorithm is taken from Qiao.Z.H's PhD thesis P41-42. 
%==========================================================================
[is_div_L,grL] = cal_surface_green_L(E,H0,H1);
[is_div_R,grR] = cal_surface_green_R(E,H0,H1);
SigmaL = HT'*grL*HT;
SigmaR = HT*grR*HT';
clear grL grR