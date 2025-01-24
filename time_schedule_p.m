function [f] = time_schedule_p(p,t,kappa)
%TIME_SCHEDULE Summary of this function goes here
%   Detailed explanation goes here

g = 1+t*(kappa^(p-1)-1);
f = kappa/(kappa-1)*(1-g^(1/(-p+1)));
 
end

