function [evecs, evals]=hspeigs(II, JJ, SS, AA, Numv)
% hou sparse eigs

W=sparse(II, JJ, SS);
Am = sparse([1:length(AA)], [1:length(AA)], AA);

nev = min(Numv, length(AA));

[evecs evals] = eigs(W, Am, nev, -1e-5);
evals = diag(evals);
%evals = evals./evals(2);

%area = sum(A);
   
%A = (1/area) * A;
   
%evals = area * evals;
   
%evecs = sqrt(area) * evecs;