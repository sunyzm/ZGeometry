function [VM, VA] = spectral_embedding(AM)

%[VM,VAM] = eigs(AM,1);
%VA = diag(VAM);

%[u v] = eig(AM);
%v = diag(v);
%[Di,IJ]=sort(v,'descend');

%VM=u(:,IJ(1));   

%VA = Di(1);

[u,s,v] = svd(AM);
VA = s(1,1);
VM = u(:,1);