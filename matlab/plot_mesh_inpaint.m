function plot_mesh_inpaint(vertex, faces, vertex_est, missing_idx)
N = size(vertex,2);
missing_percent = round(length(missing_idx) * 100 / N); 

subplot(1,2,1);
options.face_vertex_color = repmat([0.53 0.7 0.93],N,1);
options.face_vertex_color(missing_idx,:) = repmat([0 0 0],length(missing_idx), 1);
plot_mesh(vertex, faces, options)
shading faceted
title(['Original mesh with ', num2str(missing_percent), '% missing vertices']);
subplot(1,2,2);
options.face_vertex_color = repmat([0.53 0.7 0.93],N,1);
plot_mesh(vertex_est, faces, options);
shading faceted
title('Inpainted mesh');