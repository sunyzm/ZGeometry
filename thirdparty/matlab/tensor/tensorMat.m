function [vX2, mX2, score] = tensorMat(feat1, feat2, tris, numbs);

nP1 = numbs(1);
nP2 = numbs(2);
nT = numbs(3);
nNN = numbs(4);
%find the nearest neighbors
[inds, dists] = annquery(feat2, feat1, nNN, 'eps', 10); % feat2 is complete seta nd feat1 is query set
tris = int32(tris);
%build the tensor
[i j k]=ind2sub([nP2,nP2,nP2],inds);  % find multi-script representations corresponding to indices 
tmp=repmat(1:nT,nNN,1);
indH = tris(:,tmp(:))'*nP2 + [k(:)-1 j(:)-1 i(:)-1];
valH = exp(-dists(:)/mean(dists(:)));
%initiatialize X
X=1/nP2*ones(nP2,nP1);

%power iteration
[X2, score]=tensorMatching(X,[],[],[],[],indH,valH);

[vX2 mX2] = max(X2);