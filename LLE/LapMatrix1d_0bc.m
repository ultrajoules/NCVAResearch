% Computes the Matrix A such that A*F = Laplacian(F)
%
% F=N and dx = discretization

function A=LapMatrix1d_0bc(N)

A=-2*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);

return