function [ idx, coeff ] = spamsOMP( X, D, L0)
% Name: spamsOMP
%
% Description: call mexOMP to solve 
%   min_{alpha} ||x-D*alpha||_2^2  s.t. ||alpha||_0 <= L
%
% Inputs:  X: double m x 1 vector (input signal)
%          D: double m x p matrix (dictionary)
%          L: maximum number of elements in decomposition
% Outputs: idx: int L x 1 vector (indices of selected atom)
%          coeff: double L x 1 vector (coefficients)
param.L = L0;
tic;
alpha = mexOMP(X,D,param);
t = toc;
fprintf('Decomposition time(s): %f\n',t);
[idx,~,coeff] = find(alpha);
end

