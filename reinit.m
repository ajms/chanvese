function [ psi_n ] = reinit( phi, dtau, epsilon, h )
% Function to reinitialize phi as a signed distance map.
%
%    psi_n = reinit( phi, dtau, epsilon, h ) returns the
%    reinitialized levelset.
% 
%    Input:
%    - phi: levelset function to reinitialize
%    - dtau: timestep for the iterations
%    - epsilon: amout of smoothing for the sign-function
%    - h: Space step size for the approximations of the
%    derivatives.
    
    psi_n = phi;
    psi_nn = zeros(size(phi));
    while norm(psi_nn - psi_n) > h
        psi_nn = psi_n;
        psi_n = psi_n - dtau*S(phi,epsilon).*G(psi_n,phi,h);
    end
end