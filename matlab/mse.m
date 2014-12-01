function err = mse(s1, s2)
% Name: mse
% Description: return the mean squared error between s1 and s2; s1 and s2
% should both be n*m array; n is the number of elements, m is the
% number of channels.

N = size(s1,1);
sdiff = s1 - s2;
sdiff = sdiff.*sdiff;
err = sum(sum(sdiff)) / N;
end
