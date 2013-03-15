function WC = wavelet_coeff(scaling, wavelet, CLR)

nl = size(wavelet, 2);

for i = 1 : nl
	W = sparse( wavelet(i).ii, wavelet(i).jj, wavelet(i).vv );
	WC(:, i) = W*CLR;
end

