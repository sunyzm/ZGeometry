function rgb = scalar2jet(vec)
rgb = zeros(length(vec),3);
vmin = min(vec);
vmax = max(vec);
dv = vmax - vmin;

for i=1:length(vec)
    r = 1; g = 1; b = 1; 
    v = vec(i);
    if v < vmin + 0.25 * dv
        r = 0;
        g = 4 * (v - vmin) / dv;
    elseif v < vmin + 0.5 * dv
        r = 0;
        b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
    elseif v < vmin + 0.75 * dv
        r = 4 * (v - vmin - 0.5 * dv) / dv;
        b = 0;
    else 
        g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
        b = 0;
    end
    rgb(i,:)=[r g b];
end