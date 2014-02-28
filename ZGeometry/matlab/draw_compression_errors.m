function [] = draw_compression_errors( errors )
%draw_compression_errors Draw a line chart of errors against compression
%ratios
%   errors: 5*n matrix
x=errors(1,:);
y1=errors(2,:);
y2=errors(3,:);
y3=errors(4,:);
y4=errors(5,:);
%set(0,'DefaultAxesLineStyleOrder',{'-ob','--sg',':+r','-.vc'});
figure
plot(x,y1,'-ob',x,y2,'--sg',x,y3,':+r','LineWidth',2,x,y4,'-.vc')
xlabel('Compression Ratio')
ylabel('Geometric + differential error')
legend('MHB, truncation', 'MHB, S-MP', 'SGW, S-OMP', 'SGW+MHB, S-OMP')

end

