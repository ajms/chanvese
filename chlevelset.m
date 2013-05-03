function [ phi_nn ] = chlevelset( phi_n, I, lambda1, lambda2, mu, nu, dt, h )
% Update the Level Set function using the Chan-Vese Levelset.
%   chlevelset(phi_n,u0,lambda1,lambda2,mu,nu,dt,h) returns the n+1'th
%   timestep of the Level Set. Using the n'th timestep of the Level Set:
%   phi_n; the weight of exterior and interior: lambda1, lambda2; the
%   weight of the length of the contour mu; the weight of the area of the 
%   interior; timestep size dt; space step size h.

% image dimensions.
[M,N] = size(I);

% zeta and xi determine the gradient of phi.
zeta = 1./sqrt(((phi_n(2:end, 2:end-1) - phi_n(1:end-1,2:end-1)).^2/h^2 + (phi_n(1:end-1,3:end)-phi_n(1:end-1,1:end-2)).^2)./(2*h)^2);
eta = 1./sqrt(((phi_n(2:end-1,2:end)-phi_n(2:end-1,1:end-1)).^2/h^2 + (phi_n(3:end,1:end-1)-phi_n(1:end-2,1:end-1)).^2)./(2*h)^2);

%zeta = (abs(zeta)>0.001).*zeta + (abs(zeta)<=0.001).*0.001;
%eta = (abs(eta)>0.001).*eta + (abs(eta)<=0.001).*0.001;

% A,B,C,D,E are the entries of the matrix for the system to solve.
diracphi = drac(phi_n(2:end-1,2:end-1),h);
A = -mu*diracphi.*zeta(2:end,1:end)./h^2;
B = -mu*diracphi.*eta(1:end,2:end)./h^2;
C = -mu*diracphi.*zeta(1:end-1,1:end)./h^2;
D = -mu*diracphi.*eta(1:end,1:end-1)./h^2;
E = mu*diracphi.*(zeta(2:end,1:end)+eta(1:end,2:end)+zeta(1:end-1,1:end)+eta(1:end,1:end-1))./h^2 + 1/dt;

% P denotes the pentadiagonal matrix to solve.
P = spdiags(transpose([[0,A(1:end-1)];[zeros(1,M),B(1:end-M)];[C(2:end),0];[D(M+1:end),zeros(1,M)];E(1:end)]),[1 M -1 -M 0],M*N,M*N);

% build result vector from area, interior, exterior terms and boundary.
[c1 c2] = avg_intensity(phi_n(2:end-1,2:end-1), I, h);

b = phi_n(2:end-1,2:end-1)/dt + diracphi.*(-nu*ones(M,N) - lambda1*(I-c1).^2 + lambda2*(I-c2).^2);

for k=0:N-1
    P(k*N+1,k*N+1) = P(k*N+1,k*N+1) + 4/3 * C(k*N+1); % Left boundary:
    P(k*N+1,k*N+2) = P(k*N+1,k*N+2) - 1/3 * C(k*N+1); % Forward approximation
    P((k+1)*N,(k+1)*N) = P((k+1)*N,(k+1)*N) + 4/3 * A((k+1)*N); % Right boundary:
    P((k+1)*N,(k+1)*N-1) = P((k+1)*N,(k+1)*N-1) - 1/3 * A((k+1)*N); % Backward approxiNation
    P(k+1,k+1) = P(k+1,k+1) + 4/3 * D((k+1)*N); % Buttom boundary:
    P(k+2,k+1) = P(k+2,k+1) - 1/3 * D((k+1)*N); % Forward approximation
    P(end-k,end-k) = P(end-k,end-k) + 4/3 * B(end-k); % Top boundary
    P(end-(k+1),end-k) = P(end-(k+1),end-k) - 1/3 * B(end-k); % Backward approximation
    if k > 0
        P(k*N,k*N+1) = 0;
        P(k*N+1,k*N) = 0;
    end
end

% Solve linear system.bicgstab
phi_nn = zeros(M+2,N+2);
phi_nn(2:end-1,2:end-1) = reshape(P\b(:),M,N);

% Add boundary to n+1'st timestep
phi_nn(1,:) = 4/3*phi_nn(2,:) - 1/3*phi_nn(3,:);
phi_nn(end,:) = 4/3*phi_nn(end-1,:) - 1/3*phi_nn(end-2,:);
phi_nn(:,1) = 4/3*phi_nn(:,2) - 1/3*phi_nn(:,3);
phi_nn(:,end) = 4/3*phi_nn(:,end-1) - 1/3*phi_nn(:,end-2);

end