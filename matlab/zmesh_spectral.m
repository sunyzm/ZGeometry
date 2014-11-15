function [ ef1, ev1, ef2, ev2] = zmesh_spectral(mat_umbrella,mat_cot,mat_weight)
%zmesh_spectral Decompose the mesh Laplacian
%   
order = size(mat_umbrella,1);
opts.issym=1;
[ef1, ev1]=eigs(mat_umbrella,order-1,'sm',opts);
[ef2, ev2]=eigs(mat_cot, mat_weight, order-2, 'sm',opts);
end

