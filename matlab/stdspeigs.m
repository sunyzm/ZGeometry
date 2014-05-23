function [evecs, evals]=stdspeigs(II, JJ, SS, order, Numv)
% hou sparse eigs

W = sparse(II, JJ, SS, order, order);
nev = min(Numv, order-1);
opts.issym = 1;

[evecs evals] = eigs(W, nev, 'sm', opts);
evals = -diag(evals);
evecs = -evecs;
%evals = evals./evals(2);

%area = sum(A);
   
%A = (1/area) * A;
   
%evals = area * evals;
   
%evecs = sqrt(area) * evecs;