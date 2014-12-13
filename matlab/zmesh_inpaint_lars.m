function result = zmesh_inpaint_lars(dict, coord, missing_idx, lambda)
opts.lambda = lambda;
result = mesh_inpaint_lasso(dict, coord, missing_idx, opts);