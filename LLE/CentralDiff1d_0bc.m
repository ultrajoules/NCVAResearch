% Computes the Matrix B such that B*F = CentralDiff(F)
%
% F=N and dx = discretization

function B=CentralDiff1d_0bc(N)

B=+diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);

return