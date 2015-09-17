function [is_div,grL] = cal_surface_green_L(E,H0,H1)
%==========================================================================
%  Calculate the surface Green's function of left lead: gr00. 
%==========================================================================
H1_dagger = H1;
clear H1    
H1 = H1_dagger';
[is_div,grL] = cal_surface_green_R(E,H0,H1);