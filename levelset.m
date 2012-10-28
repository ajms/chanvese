function [ phi_nn ] = levelset( phi_n, u0, lambda1, lambda2, mu, nu, dt, h, doreinit )
%LEVELSET Summary of this function goes here
%   Detailed explanation goes here
    
    [c1 c2] = avg_intensity(phi_n, u0, h);
    
    Dxphi = (phi_n(3:end, 2:end-1) - phi_n(1:end-2, 2:end-1))/(2*h);
    Dyphi = (phi_n(2:end-1, 3:end) - phi_n(2:end-1, 1:end-2))/(2*h);
    Dxxphi = (phi_n(3:end, 2:end-1) - 2*phi_n(2:end-1, 2:end-1) + phi_n(1:end-2, 2:end-1))/(h^2);
    Dyyphi = (phi_n(2:end-1, 3:end) - 2*phi_n(2:end-1, 2:end-1) + phi_n(2:end-1, 1:end-2))/(h^2);
    Dxyphi = (phi_n(3:end, 3:end) - phi_n(3:end, 1:end-2) - phi_n(1:end-2, 3:end) + phi_n(1:end-2, 1:end-2))/(4*h^2);
    
    g = sqrt(Dxphi.^2 + Dyphi.^2);
    
    if doreinit == 0
        g(g<0.5) = 1;
    end
    
    k = (Dxphi.^2.*Dyyphi + Dyphi.^2.*Dxxphi - 2*Dxphi.*Dyphi.*Dxyphi)./g.^3;
    
    if doreinit == 0
        k(k > 1/h) = 1/h;
        k(k < -1/h) = -1/h;
    end
    
    phi_nL = phi_n(2:end-1, 2:end-1);
    u0L = u0(2:end-1, 2:end-1);
    
    phi_nn = phi_n;
    
    % PDE step
    phi_nn(2:end-1, 2:end-1) = drac(phi_nL, h).*(mu*k - nu + lambda2*(u0L - c2).^2 - lambda1*(u0L - c1).^2)*dt + phi_nL;
    
    % Apply boundary conditions
    phi_nn(2:end-1, 1) = phi_nn(2:end-1, 2);
    phi_nn(2:end-1, end) = phi_nn(2:end-1, end-1);
    phi_nn(1, 2:end-1) = phi_nn(2, 2:end-1);
    phi_nn(end, 2:end-1) = phi_nn(end-1, 2:end-1);
    
    % Re-initialization
    if doreinit == 1
        phi_nn = reinit(phi_nn);
    end
end
