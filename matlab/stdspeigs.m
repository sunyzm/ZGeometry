function [evecs, evals]=stdspeigs(II, JJ, SS, order, Numv)
% sparse standard eigendecomposition

Ls = sparse(II, JJ, SS, order, order);
nev = min(Numv, order-1);
opts.issym = 1;

% solve
[evecs evals] = eigs(Ls, nev, 'sm', opts);
evals = -diag(evals);
evecs = -evecs;
