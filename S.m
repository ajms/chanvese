function [ sign ] = S( phi, epsilon )
% Calculates a smoothed sign function
%   sign = S(phi,epsilon) returns a smooth sign function
%
%   Input:
%   - phi: levelset function
%   - epsilon: amount of smoothing.
    
    sign = phi./sqrt(phi.^2+epsilon^2);
end