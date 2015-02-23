function [result,alpha, err] = zmesh_inpaint_lars(dict, coords, missing_idx, lambda)
tic

N = size(coords, 1);
projmask = true(N,3);
projmask(missing_idx,:) = false;

param.mode = 1;
param.lambda = lambda;
alpha = mexLassoMask(coords, dict, projmask, param);
result = dict * alpha;

elapsed = toc;
err = mse(coords(missing_idx,:), result(missing_idx,:));
missing_percent = round(length(missing_idx) * 100 / N); 
disp(['Mesh inpainting error (LARS, ', num2str(missing_percent), '%): ', num2str(err)]);
disp(['Elapsed time (LARS, ', num2str(missing_percent), '): ', num2str(elapsed)]);
 
end