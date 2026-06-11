function [c,cx,cy]=IceVelocity2CalvingRate(CtrlVar,MUA,F,LSF,u,v,nx,ny)

% Calculates the calving rate that equals the normal velocity to the calving front
%
%   $$ c = v \cdot n $$
%

narginchk(6,8)
nargoutchk(1,3)

if nargin<7
    [dLSFdx,dLSFdy]=calcFEderivativesMUA(LSF,MUA,CtrlVar);
    [nx,ny]=ProjectFintOntoNodes(MUA,dLSFdx,dLSFdy);
    Norm=vecnorm([nx ny],2,2);
    nx=nx./(Norm+eps) ;
    ny=ny./(Norm+eps) ;
end

if isempty(u) || isempty(v)
    c=[]; cx=[] ; cy=[] ; 
    return
end

c=-(u.*nx+v.*ny);  % the minus is because the normal of LSF points upstream

%
%
%



if nargout>1
    cx=c.*nx ; cy=c.*ny;

end


end