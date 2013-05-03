function [ rephi ] = reinit( phi, dtau, epsilon, h )
    psi_n = phi;
    psi_nn = zeros(size(phi));
    while norm(psi_nn - psi_n) > h
        fprintf('Reiniterror: %f\n', norm(psi_nn-psi_n));
        psi_nn = psi_n;
        psi_n = psi_n - dtau*S(phi,epsilon).*G(psi_n,h);
        surf(psi_nn);
        pause(0.1);
    end
end