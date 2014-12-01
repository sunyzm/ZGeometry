function [ef, ev] = zmesh_spectral_std(varargin)
%zmesh_spectral_std Decompose the standard mesh Laplacian
%   
mat_sym = varargin{1};
order = size(mat_sym,1);
opts.issym = 1;
if nargin == 1
    [ef, ev] = eigs(mat_sym, order-1, 'sm', opts);
else
    nev = varargin{2};
    [ef, ev] = eigs(mat_sym, min(order-1,nev), 'sm', opts);
end
end

