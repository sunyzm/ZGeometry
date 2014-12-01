function s = random_subset(N, n)
% Name: random_subset
% Descrption: return a random n indices from 1:N

A = 1:N;
s = zeros(1,n);
for i = 1:n
    len = length(A);
    idx = randi(len);
    s(i) = A(idx);
    A = horzcat(A(1:idx-1),A(idx+1:len));
end
s = sort(s);
