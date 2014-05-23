function [evecs, evals]=hspeigs(II, JJ, SS, AA, Numv)
% hou sparse eigs

order = length(AA);
W = sparse(II, JJ, SS, order, order);
Am = sparse(diag(AA));
nev = min(Numv, order-1);
opts.issym = 1;

[evecs evals] = eigs(W, Am, nev, 'sm', opts);
evals = -diag(evals);
evecs = -evecs;
%evals = evals./evals(2);

%area = sum(A);
   
%A = (1/area) * A;
   
%evals = area * evals;
   
%evecs = sqrt(area) * evecs;