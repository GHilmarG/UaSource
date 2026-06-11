function [daMin,ia,dbMin,ib,D2]=DistanceCloudAtoCloudB(xa,ya,xb,yb)
% [daMin,ia,dbMin,ib,D2]=DistanceCloudAtoCloudB(xa,ya,xb,yb)
% finds distance between all points in point cloud A to all points in cloud B
% also the closest point in one cloud to a given point in the other.
%
% d     : a distance matrix
% daMin : a vector giving the min distance of every point in cloud A to a cloud B
% ia    : a vector giving the number of the point in cloud B closest to the points of cloud A
%
% Example: daMin(1) will give the minimum distance from [xa(1) ya(1)] to the points in cloud B
%          and ia(1) will be the number of that node in cloud B
%
% SEE ALSO:
% RangeSearch
%
%
%
xa=xa(:) ;  ya=ya(:) ;  xb=xb(:) ;  yb=yb(:) ;

[Xa,Xb]=ndgrid(xa,xb) ;
[Ya,Yb]=ndgrid(ya,yb);


D2=(Xa-Xb).^2+(Ya-Yb).^2;  % this is symmetricc so I am calculating a bit too much here

%D2= bsxfun(@minus,Xa,Xb).^2 + bsxfun(@minus,Ya,Yb).^2; % slower! (2016b)


[daMin,ia]=min(D2,[],2);
daMin=sqrt(daMin(:)); ia=ia(:);

if nargout>2
    [dbMin,ib]=min(D2,[],1);
    dbMin=sqrt(dbMin(:)); ib=ib(:);
else
    dbMin=[]; D2=[];
end

end

