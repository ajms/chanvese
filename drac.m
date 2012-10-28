function [ d ] = drac( z, eps )
%DIRAC Summary of this function goes here
%   Detailed explanation goes here

    d = eps./(pi.*(z.^2 + eps^2));
end

