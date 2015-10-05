
addpath('ann_mwrapper');
addpath('mex');
%number of point in image 1
nP1=50;
%number of point in image 2
nP2=nP1;

%randomly generate them
P1=randn(2,nP1);

%generate modified version of points 1
scale= 0.5+rand();
theta = 0.5*rand();
Mrot = [cos(theta) -sin(theta) ; sin(theta) cos(theta) ];
P2=Mrot*P1*scale+0.05*randn(size(P1));


%number of used triangles (results can be bad if too low)
nT=nP1*50;
t1=floor(rand(3,nT)*nP1);
while 1
  probFound=false;
  for i=1:3
    ind=(t1(i,:)==t1(1+mod(i,3),:));
    if(nnz(ind)~=0)
      t1(i,ind)=floor(rand(1,nnz(ind))*nP1);
      probFound=true;
    end
  end
  if(~probFound)
    break;
  end
end

%generate features
t1=int32(t1);
[feat1,feat2] = mexComputeFeature(P1,P2,int32(t1),'simple');


%number of nearest neighbors used for each triangle (results can be bad if
%too low)
nNN=300;

%find the nearest neighbors
[inds, dists] = annquery(feat2, feat1, nNN, 'eps', 10);

%build the tensor
[i j k]=ind2sub([nP2,nP2,nP2],inds);
tmp=repmat(1:nT,nNN,1);
indH = t1(:,tmp(:))'*nP2 + [k(:)-1 j(:)-1 i(:)-1];
valH = exp(-dists(:)/mean(dists(:)));
%initiatialize X
X=1/nP2*ones(nP2,nP1);

%power iteration
[X2, score]=tensorMatching(X,[],[],[],[],indH,valH);



%draw
if 1
  figure(1);
  imagesc(X2);
  figure(2);
  hold on;
  plot(P1(1,:),P1(2,:),'r x');
  plot(P2(1,:),P2(2,:),'b o');
  [tmp match] = max(X2);
  for p=1:nP1
    plot([P1(1,p),P2(1,match(p))],[P1(2,p),P2(2,match(p))],'k- ');
  end
end










