function hks = prolongHeatKernel(proMat, coarseHK)
[m,]=size(proMat);
hks = zeros(m,1);
for k = 1:m
    %rowi = proMat(k,:);
    hks(k) = proMat(k,:) * coarseHK * proMat(k,:)';
end
