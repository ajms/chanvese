function [ dice ] = dice( phi, truth )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    seg1 = phi;
    seg2 = truth;
    seg1(seg1 > 0) = 1;
    seg1(seg1 <= 0) = -1;
    seg2(seg2 < 100) = 1;
    seg2(seg2 >= 100) = 0;
    dice = 2*(sum(sum(seg1 == seg2)))/(sum(sum(seg1 == 1)) + sum(sum(seg2)));
end

