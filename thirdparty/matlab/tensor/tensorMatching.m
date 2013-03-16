function [X2 score] = tensorMatching(X,indH1,valH1,indH2,valH2,indH3,valH3,nIter,sparsity,stoc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [X2 score] = tensorMatching(X,indH1,valH1,indH2,valH2,indH3,valH3[,nIter,sparsity,stoc])
%input:
%  X / N1 x N2 matrix / initial assignment matrix
%  indHi / Nt x i matrix / indices in sparse tensor of order i
%  valHi / Nt x 1 / values in sprase tensor of order i
%  nIter / scalar integer / number of iteration for the power method
%  sparsity / boolean / should the algorithm use the sparsity trick as
%                         discribed in the article
%  stoc / boolean / 0 for no contraint on X, 1 for stochastic
%                         constraint, 2 for doubly stochastic
%  
%output:
%  X2 / N1 x N2 matrix / output assignment matrix
%  score / scalar / score of the matching (high means good matching)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(indH1) && isempty(valH1)
  indH1=int32(zeros(0,1));
  valH1=zeros(0,1);
end
if isempty(indH2) && isempty(valH2)
  indH2=int32(zeros(0,2));
  valH2=zeros(0,1);
end
if isempty(indH3) && isempty(valH3)
  indH3=int32(zeros(0,3));
  valH3=zeros(0,1);
end

if nargin<8
  nIter=100;
end
if nargin<9
  sparsity=1;
end
if nargin<10
  stoc=2;
end
[X2, score]=mexTensorMatching(double(X),indH1,valH1,indH2,valH2,indH3,valH3,nIter,sparsity,stoc);

