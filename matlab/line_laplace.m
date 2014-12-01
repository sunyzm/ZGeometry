function L = line_laplace(n)
L=spalloc(n,n,3*n);
L(1,1) = 1; L(1,2) = -1;
L(n,n) = 1; L(n,n-1) = -1;
for i = 2:n-1
    L(i,i) = 2;
    L(i,i-1) = -1;
    L(i,i+1) = -1;
end
end
