function [ phi_nn ] = chlevelset( phi_n, u0, lambda1, lambda2, mu, nu, dt, h )
% Update the Level Set function using the Chan-Vese Levelset.
%   chlevelset(phi_n,u0,lambda1,lambda2,mu,nu,dt,h) returns the n+1'th
%   timestep of the Level Set. Using the n'th timestep of the Level Set:
%   phi_n; the weight of exterior and interior: lambda1, lambda2; the
%   weight of the length of the contour mu; the weight of the area of the 
%   interior; timestep size dt; space step size h.

% image dimensions.
m = size(phi_n,1)-2;

% zeta and xi determine the gradient of phi.
zeta = sqrt(((phi_n(3:end, 2:end-1) - phi_n(2:end-1,2:end-1)).^2 + (phi_n(2:end-1,3:end)-phi_n(2:end-1,1:end-2)).^2)/h^2);

xi = sqrt(((phi_n(2:end-1,3:end)-phi_n(2:end-1,2:end-1)).^2+(phi_n(3:end,2:end-1)-phi_n(1:end-2,2:end-1)).^2)/h^2);


% check ./ vs. /
% C,D,E are the entries of the matrix for the system to solve.
C = reshape(drac(phi_n(2:end-1,2:end-1),h)*mu./(h^2*xi),1,m^2); % y
D = reshape(drac(phi_n(2:end-1,2:end-1),h)*mu./(h^2*zeta),1,m^2); % x
E = reshape(drac(phi_n(2:end-1,2:end-1),h)*(-2*mu)./(h^2*(xi-zeta))-1/dt,1,m^2);

fprintf('%d\n',size(C,2),size(E,2),size(D,2))

% proper formulation for sparse matrices.
% A denotes the matrix to solve.
A = sparse(diag(E(1:end),0)) + sparse(diag(D(1:end-1),1)) + sparse(diag(D(2:end),-1)) + sparse(diag(C(m+1:end),-m)) + sparse(diag(C(1:end-m),m));

% missing proper boundaries
% Boundary conditions.
gl = 1;
gr = 1;
gt = 1;
gb = 1;

% Build result vector from area, interior, exterior terms and boundary.
u0L = u0(2:end-1,2:end-1);
[c1 c2] = avg_intensity(phi_n, u0, h);

b = reshape(phi_n(2:end-1,2:end-1),m^2,1)/dt + nu*ones(m^2,1) + lambda1*(reshape(u0L,m^2,1)-c1).^2 - lambda2*(reshape(u0L,m^2,1)-c2).^2;

for j=1:m
    b(j) = b(j) + C(j)*gb;
    b(end-j) = b(end-j) + C(end-j)*gt;
    b(j*m) = b(j*m) + D(j*m)*gr;
    b((j-1)*m+1) = b((j-1)*m+1) + D((j-1)*m+1)*gl;
    if j > 1
        A((j-1)*m,(j-1)*(m+1)) = 0;
        A((j-1)*(m+1),(j-1)*m) = 0;
    end
end

% eventually find another numerical method to solve the system.
% Solve matrix system using conjugent gradient method (pcg).
phi_nL = linsolve(full(A),b);

phi_nn = ones(m+2,m+2);
phi_nn(2:end-1,2:end-1) = reshape(phi_nL,m,m);
phi_nn(1,2:end-1) = phi_nn(2,2:end-1);
phi_nn(end,2:end-1) = phi_nn(end-1,2:end-1);
phi_nn(2:end-1,1) = phi_nn(2:end-1,2);
phi_nn(2:end-1,end) = phi_nn(2:end-1,end-1);

end