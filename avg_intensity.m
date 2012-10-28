function [ c1, c2 ] = avg_intensity ( phi, u0, h )
%AVERAGE Summary of this function goes here
%   Detailed explanation goes here
    c1 = sum(sum(u0.*hside(phi, h)))/sum(sum(hside(phi, h)));
    c2 = sum(sum(u0.*(1 - hside(phi, h))))/sum(sum(1 - hside(phi, h)));
end

