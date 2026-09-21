

function [xeff,dxeffdx,d2xeffdx2]=SmoothFloor(x,xmin,w)

% Smooth replacement for xeff=max(x,xmin).
%   xeff > xmin strictly, for any real x
%   xeff -> x + w^2/(4(x-xmin))  for x >> xmin

if w<=0 
    error("SmoothFloor:zeroWidth","Transition width must be strictly positive.") ; 
end

d = x - xmin ;
r = sqrt(d.*d + w*w) ;
xeff = xmin + (d+r)/2 ;
if nargout>1 
    dxeffdx   = (1 + d./r)/2 ; 
end

if nargout>2 
    d2xeffdx2 = w*w./(2*r.^3) ;
end

end