function [ idx, coeff, nnz ] = spamsLasso( X, D, L1)
% Description: call mexLASSO to solve 
%   min_{alpha} 0.5||x-Dalpha||_2^2 + lambda||alpha||_1
%
% Inputs:  X: double m x 1 vector (input signal)
%          D: double m x p matrix (dictionary)
%          L1: parameter
% Outputs: idx: int L x 1 vector (indices of selected atom)
%          coeff: double L x 1 vector (coefficients)
param.lambda=L1;
param.mode=2;
alpha=mexLasso(X,D,param);
[idx, col, coeff]=find(alpha);
nnz=length(idx);
end