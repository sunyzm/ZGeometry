function [] = draw_compression_errors( errors )
%draw_compression_errors Draw a line chart of errors against compression
%ratios
%   
x=errors(1,:);
y1=errors(2,:);
y2=errors(3,:);
y3=errors(4,:);
y4=errors(5,:);

figure
hold on
plot(x,y1,'Color','k','LineStyle','-','Marker','o','LineWidth',1)
plot(x,y2,'Color','r','LineStyle','--','Marker','s','LineWidth',1)
plot(x,y3,'Color',[0,0.5,0],'LineStyle',':','Marker','+','LineWidth',2)
plot(x,y4,'Color','b','LineStyle','-.','Marker','x','LineWidth',1.5)
xlim([0.03,0.82])
xlabel('Compression Ratio','FontSize',13)
ylabel('Geometric + Differential Error','FontSize',13)
hleg = legend('MHB, truncation', 'MHB, S-MP', 'SGW, S-OMP', 'SGW+MHB, S-OMP')
set(hleg,'FontSize',12,'Location','NorthEast')
grid on
hold off
end

