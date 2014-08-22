function evals = sparse_bi_conj_grad(II,JJ,SS,m,n,vy)
%m,n: size of sparse lefthand matrix
%II,JJ,SS: define the sparse matrix m*n
%vy: righthand vector m*1
%evals: return value n*1
A = sparse(II, JJ, SS, m, n);
evals = cgs(A'*A, A'*vy); 