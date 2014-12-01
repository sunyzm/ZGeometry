function [ef, ev] = zmesh_spectral_gen(varargin)
%zmesh_spectral Decompose the mesh Laplacian
%   
mat_sym = varargin{1};
mat_weight = varargin{2};
order = size(mat_sym,1);
opts.issym = 1;
if nargin == 2
    [ef, ev] = eigs(mat_sym, mat_weight, order-2, 'sm', opts);
else
    nev = varargin{3};
    [ef, ev] = eigs(mat_sym, mat_weight, min(nev, order-2), 'sm', opts);
end
end

