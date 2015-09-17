function f = cal_fermi(E,mu,Ef,kT)
%==========================================================================
%   calculate Fermi distribution functions. Round to zero if
%   :  (E+mu-Ef)/kT > 40.0
%==========================================================================

if (E+mu-Ef)/kT > 40.0
    f = 0.0;
else
    f = 1.0/(exp((E+mu-Ef)/kT) + 1);
end
