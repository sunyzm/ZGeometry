function h = zmesh_draw(verts, faces, colors)

if size(colors,1) == size(verts,1)
    shading_type = 'interp';
else
    shading_type = 'flat';
end

h = patch('vertices',verts,'faces',faces, 'FaceVertexCData', colors, 'FaceColor', shading_type);
lighting phong;
camproj('perspective');
axis off;
axis tight;
axis equal;
shading interp;
camlight;

end
