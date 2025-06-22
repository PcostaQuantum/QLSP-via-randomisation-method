function [s] = func_s(j,q,kappa)
%FUNC_S Summary of this function goes here
%
% Input: 
%   j: index of the current step
%   q: overall number of steps
%   kappa: condition number

va = sqrt(2)*kappa / sqrt(1+kappa^2) * log(kappa*sqrt(1+kappa^2) - kappa^2); 
vb = sqrt(2)*kappa / sqrt(1+kappa^2) * log(sqrt(1+kappa^2) + 1); 

v = va + j*(vb-va)/q; 

s = (-kappa^2 * exp(-v*sqrt(kappa^2+1)/sqrt(2)/kappa) ...
    + exp(v*sqrt(kappa^2+1)/sqrt(2)/kappa) + 2*kappa^2) ...
    / 2 / (kappa^2+1); 


end

