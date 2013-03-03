function [ k,x ] = jacobi( A, b, M )
%jacobi(A,b) solves the linear equation system with the iterative
%Jacobi method
%   Input arguments: nxn matrix A and nx1 vector b.

n = size(b,1);
x = ones(n,1);
u = zeros(n,1);
error = 1;
k = 0;
while error > 10^(-6) && k < M
    for i=1:n
        j = 1:n;
        j(i) = [];
        u(i) = (b(i)-sum(A(i,j)*x(j)))/A(i,i);
    end
    error = norm(u-x,2);
    x = u;
    k = k+1;
end
end

