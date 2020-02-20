function [diffDu,diffDh]=CalcIncrementsNorm(CtrlVar,MUA,L,Inodes,u,du,v,dv,h,dh)

narginchk(8,10)
nargoutchk(1,2)


du(Inodes)=0;
dv(Inodes)=0;

u=abs(u)+CtrlVar.SpeedZero ; 
v=abs(v)+CtrlVar.SpeedZero ; 
du=du./u ;
dv=dv./v; 


% This will be equal to x if all elements of dub and dvb are equal to x
% It is an area integral over the residuals and then divided by the area
% The units are the same as those of dub, etc
diffDu=sqrt((du'*MUA.M*du+dv'*MUA.M*dv)/2/MUA.Area);

if nargout>1
    dh(Inodes)=0;
    h=h+CtrlVar.ThickMin; 
    dh=dh./h ; 
    diffDh=sqrt( (dh'*MUA.M*dh)/MUA.Area);
else
    diffDh=[];
end

end
%%