function [ absgrad ] = G( psi, phi, h )
% Calculates the gradient of psi-1.
%   absgrad = G(psi,phi,h) returns the absolute gradient -1.
%   
%   Input: 
%   - psi: Levelset in the n'th iteration in the reinitialization
%   - phi: Not reinitialized levelset function
%   - h: space step size
    
    a = (psi(2:end-1,2:end-1)-psi(1:end-2,2:end-1))/h;
    b = (psi(3:end,2:end-1)-psi(2:end-1,2:end-1))/h;
    c = (psi(2:end-1,2:end-1)-psi(2:end-1,1:end-2))/h;
    d = (psi(2:end-1,3:end)-psi(2:end-1,2:end-1))/h;
    absgrad = zeros(size(psi));
    for i=1:size(phi,1)-2
        for j=1:size(phi,2)-2
            if phi(i+1,j+1) > 0
                absgrad(i+1,j+1) = sqrt(max((a(i,j)>0)*a(i,j)^2,(b(i,j)<0)*b(i,j)^2))+sqrt(max((c(i,j)>0)*c(i,j)^2,(d(i,j)<0)*d(i,j)^2))-1;
            elseif phi(i+1,j+1) < 0
                absgrad(i+1,j+1) = sqrt(max((a(i,j)<0)*a(i,j)^2,(b(i,j)>0)*b(i,j)^2))+sqrt(max((c(i,j)<0)*c(i,j)^2,(d(i,j)>0)*d(i,j)^2))-1;
            end
        end
    end
end
