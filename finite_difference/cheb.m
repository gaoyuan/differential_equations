function [D,x]=cheb(N)
% generating Chebychev Pts and its differential matrix
if N==0, d=0; x=1; return; end
x = cos(pi*(0:N)/N)';
c = [2;ones(N-1,1);2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X - X';
D = (c*(1./c)')./(dX+eye(N+1)); % off diagonals
D = D - diag(sum(D')); % diagonals
