function [ c1, c2 ] = avg_intensity ( phi, I, h )
% Function to compute the average intensities of the segments
%
%   [c1, c2]  = avg_intensity ( phi, u0, h ) returns the average
%   intensity of the inner segment (c1) and the outer segment (c2)
%
%   Input:
%   - phi: levelset function
%   - I: image
%   - h: space step size
    
    c1 = sum(sum(I.*hside(phi, h)))/sum(sum(hside(phi, h)));
    c2 = sum(sum(I.*(1 - hside(phi, h))))/sum(sum(1 - hside(phi, h)));
end

