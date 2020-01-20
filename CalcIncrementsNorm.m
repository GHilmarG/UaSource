function [diffDu,diffDh]=CalcIncrementsNorm(CtrlVar,MUA,L,Inodes,dub,dvb,dh)

dub(Inodes)=0;
dvb(Inodes)=0;
dh(Inodes)=0; 


% This will be equal to x if all elements of dub and dvb are equal to x
% It is an area integral over the residuals and then divided by the area
% The units are the same as those of dub, etc
diffDu=sqrt((dub'*MUA.M*dub+dvb'*MUA.M*dvb)/2/MUA.Area);
diffDh=sqrt( (dh'*MUA.M*dh)/MUA.Area);

end
%%