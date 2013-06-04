function [ H ] = hside( z, eps )
% Calculates a smoothed heavyside function
%   H = hside(z,eps) returns a smooth heavyside function
%
%   Input: 
%   - z: vector/matrix of z-values
%   - eps: amount of smoothing

    H = 0.5*(1 + 2/pi*atan(z./eps));
end

