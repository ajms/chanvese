function [ k,x ] = jacobi( A, b, M )
%jacobi(A,b) solves the linear equation system with the iterative
%Jacobi method
%   Input arguments: nxn matrix A and nx1 vector b.

n = size(b,1);
x = ones(n,1);
u = zeros(n,1);
error = 1;
k = 0;

for i = 1:n
    j = 1:n;
    j(i) = [];
    B = abs(A(i,j));
    Check(i) = abs(A(i,i)) - sum(B); % Is the diagonal value greater than the remaining row values combined?
    if Check(i) < 0
        fprintf('The matrix is not strictly diagonally dominant at row %2i\n\n',i)
    end
end

while error > 10^(-2) && k < M
    Ax = bsxfun(@times,A,reshape(x,1,n));
    u = (b-sum(Ax-diag(diag(Ax),0),2))./diag(A);
    error = norm(u-x,2);
    x = u;
    k = k+1;
end
end

