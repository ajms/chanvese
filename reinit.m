function [ re_phi ] = reinit( phi )
%REINIT Summary of this function goes here
%   Detailed explanation goes here

    re_phi = double((phi > 0).*(bwdist(phi < 0) - 0.5) - (phi < 0).*(bwdist(phi > 0) - 0.5));
end

