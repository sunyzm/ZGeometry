function [evecs, evals]=hspeigs(II, JJ, SS, AA, Numv)
% sparse general eigendecomposition

order = length(AA);
Ls = sparse(II, JJ, SS, order, order);
W = sparse(diag(AA));
nev = min(Numv, order-2);
opts.issym = 1;

% solve
[evecs evals] = eigs(Ls, W, nev, 'sm', opts);
evals = -diag(evals);
evecs = -evecs;