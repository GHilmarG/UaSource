function [d,in,on]=SignedDistance(p,Boundary,ds)

%%
% Calculates signed distance from boundary.
%
% [d,in,on]=SignedDistance(points,Boundary,ds)
%
% Inputs:
% p        : Points to be tested as an Nx2 array [x1 x2 ; yx2 y2 ; etc]
% Boundary : Boundary points as an Nx2 array [x1 x2 ; yx2 y2 ; etc]
% ds       : (optional) if given, then distance will only be calculated for 
%             the subset of points in p that are within a square given by
%             xmin-ds to xmax+ds and ymin-ds to ymax+ds
%
% Outputs:
% d         : signed distance
% in        : true for points inside Boundary
% on        : true for points on Boundary
%
%

N =size(p,1);
d=zeros(N,1)+NaN;

if nargin>2
    xmin=min(Boundary(:,1)) ; xmax=max(Boundary(:,1));
    ymin=min(Boundary(:,2)) ; ymax=max(Boundary(:,2));
    I=p(:,1)> (xmin-ds) & p(:,1) <(xmax+ds) ...
        & p(:,2)> (ymin-ds) & p(:,2) <(ymax+ds) ;
    I=find(I);
else
    I=1:N;
end


% first calculate distance
for k=1:numel(I)
    temp=((p(I(k),1)-Boundary(:,1)).^2+ (p(I(k),2)-Boundary(:,2)).^2);
    d(I(k))=sqrt(min(temp));
end

% then give it the right sign
[in,on] = inpoly(p,Boundary); % much faster than matlab inpolygon!
d(~in)=-d(~in);



return

end
