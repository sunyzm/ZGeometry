function evals=cgls(A,b)

evals=cgs(A'*A, A'*b)