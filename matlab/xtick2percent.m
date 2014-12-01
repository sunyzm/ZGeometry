% Convert y-axis values to percentage values by multiplication
     a = cellstr(num2str(get(gca,'xtick')'*100)); 
% Create a vector of '%' signs
     pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
     new_xticks = [char(a),pct];
% 'Reflect the changes on the plot
     set(gca,'xticklabel',new_xticks) 