function [ H ] = hside( z, eps )
%HEAVISIDE Summary of this function goes here
%   Detailed explanation goes here

    H = 0.5*(1 + 2/pi*atan(z./eps));
end

