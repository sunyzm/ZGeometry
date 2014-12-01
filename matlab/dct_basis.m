function basis = dct_basis(N)

basis = zeros(N,N);
basis(:,1) = sqrt(1./N) * ones(N,1);

dctf=@(k,n) sqrt(2./N) * cos(pi*(k-1)*(2*n-1) / (2*N));

for i = 2:N
    basis(:,i) = dctf(i,1:N);
end
