function [ d ] = drac( z, eps )
% Calculates a smoothed dirac function
%   d = drac(z,eps) returns a smooth dirac function
%
%   Input: 
%   - z: vector/matrix of z-values
%   - eps: amount of smoothing
    
    d = eps./(pi.*(z.^2 + eps^2));
end

