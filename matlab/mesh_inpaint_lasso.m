function result = mesh_inpaint_lasso(dict, coord, missing_idx, opts)

N = size(coord, 1);
missing_percent = round(length(missing_idx) * 100 / N); 
projmask = true(N,3);
projmask(missing_idx,:) = false;

param.mode = 1;
param.lambda = opts.lambda;
alpha = mexLassoMask(coord, dict, projmask, param);
result = dict * alpha;

err = mse(coord(missing_idx,:), result(missing_idx,:));
disp(['Mesh inpainting error (lasso, ', num2str(missing_percent), '%): ', num2str(err)]);