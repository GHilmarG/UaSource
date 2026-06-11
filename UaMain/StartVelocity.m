function F=StartVelocity(CtrlVar,MUA,BCs,F)

% [F.ub,F.vb,F.ud,F.vd]=StartVelocity(CtrlVar,MUA,BCs,F)

F.ub=zeros(MUA.Nnodes,1) ;
F.vb=zeros(MUA.Nnodes,1) ;
F.ud=zeros(MUA.Nnodes,1) ;
F.vd=zeros(MUA.Nnodes,1) ;




F.ub(BCs.ubFixedNode)=BCs.ubFixedValue;
F.vb(BCs.vbFixedNode)=BCs.vbFixedValue;
    
F.ud(BCs.udFixedNode)=BCs.udFixedValue;
F.vd(BCs.vdFixedNode)=BCs.vdFixedValue;

end