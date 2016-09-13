function [ub,vb,ud,vd]=StartVelocity(CtrlVar,MUA,BCs,ub,vb,ud,vd,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m)

ub=zeros(MUA.Nnodes,1) ;
vb=zeros(MUA.Nnodes,1) ;
ud=zeros(MUA.Nnodes,1) ;
vd=zeros(MUA.Nnodes,1) ;




ub(BCs.ubFixedNode)=BCs.ubFixedValue;
vb(BCs.vbFixedNode)=BCs.vbFixedValue;
    
ud(BCs.udFixedNode)=BCs.udFixedValue;
vd(BCs.vdFixedNode)=BCs.vdFixedValue;

end