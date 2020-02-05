function [diffDu,diffDh]=CalcIncrementsNorm(CtrlVar,MUA,L,Inodes,du,dv,dh)

narginchk(6,7)
nargoutchk(1,2)


du(Inodes)=0;
dv(Inodes)=0;



% This will be equal to x if all elements of dub and dvb are equal to x
% It is an area integral over the residuals and then divided by the area
% The units are the same as those of dub, etc
diffDu=sqrt((du'*MUA.M*du+dv'*MUA.M*dv)/2/MUA.Area);

if nargout>1
    dh(Inodes)=0;
    diffDh=sqrt( (dh'*MUA.M*dh)/MUA.Area);
else
    diffDh=[];
end

end
%%