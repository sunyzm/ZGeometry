function WC = wavelet_diffusion(scaling, wavelet, CLR)

nl = size(wavelet, 2);

for i = 1 : nl
	W = sparse( scaling(i).ii, scaling(i).jj, scaling(i).vv );
	WC(:, i) = W*CLR;
end

