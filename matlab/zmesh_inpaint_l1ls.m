function [ result, alpha, err ] = zmesh_inpaint_l1ls( dict, coords, missing_idx, lambda, rel_tol)
tic

N = size(coords, 1);
remaining_idx = setdiff(1:N, missing_idx);
A = dict(remaining_idx,:);
y0 = coords(remaining_idx, :);
alpha = zeros(size(dict,2), 3);

quite = true;
[alpha(:,1), ~] = l1_ls(A, y0(:,1), lambda, rel_tol, quite);
[alpha(:,2), ~] = l1_ls(A, y0(:,2), lambda, rel_tol, quite);
[alpha(:,3), ~] = l1_ls(A, y0(:,3), lambda, rel_tol, quite);
result = dict * alpha;

elapsed = toc;
err = mse(coords(missing_idx,:), result(missing_idx,:));
missing_percent = round(length(missing_idx) * 100 / N); 
disp(['lambda = ', num2str(lambda), ' ,Mesh inpainting error (l1_ls, ', num2str(missing_percent), '%): ', num2str(err)]);
disp(['Elapsed time (l1_ls, ', num2str(missing_percent), '): ', num2str(elapsed)]);

end

