function [Cov,iCov]=SparseCovarianceDistanceMatrix(x,y,Err,Sigma,DistanceCutoff)

%%
%save TestSave x y Err Sigma DistanceCutoff

%      i = 1:8;
% matrix = abs(bsxfun(@minus,i',i));
% covariance = repmat(.5,8,8).^matrix;
%
n=length(x);
N=10*n;  % only a rough initial guess
j=zeros(N,1)+NaN ; i=zeros(N,1)+NaN ; v=zeros(N,1)+NaN;
k=0;
for I=1:n-1 % upper-right
    
    dist=sqrt((x(I)-x(I+1:end)).^2+(y(I)-y(I+1:end)).^2);
    ind=find(dist<DistanceCutoff);
    l=length(ind);
    if k+l>N
        j=[j;zeros(N,1)]; i=[i;zeros(N,1)]; v=[v;zeros(N,1)];
        N=length(j);
    end
    j(k+1:k+l)=ind+I;
    i(k+1:k+l)=ind*0+I;
    v(k+1:k+l)=dist(ind);
    k=k+l;
end

i(k+1:k+n)=1:n ; j(k+1:k+n)=1:n ; v(k+1:k+n)=0 ; k=k+n; % main diagonal
v=exp(-v/Sigma)*Err^2;



Cov=sparse(i(1:k),j(1:k),v(1:k),n,n);
Cov=(Cov+Cov')/2; % Cov(1:n,1:n)=Cov(1:n,1:n)/2;



%%

if nargout==2
    iCov=inv(Cov);
end

end


