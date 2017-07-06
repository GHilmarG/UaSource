function [p,plb,pub]=InvValues2p(CtrlVar,InvValues)







if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
    
    pA=log10(InvValues.AGlen);
    lbA=log10(CtrlVar.AGlenmin)+zeros(size(pA));
    ubA=log10(CtrlVar.AGlenmax)+zeros(size(pA));
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
    
    pA=InvValues.AGlen;
    lbA=CtrlVar.AGlenmin+zeros(size(pA));
    ubA=CtrlVar.AGlenmax+zeros(size(pA));
    
elseif ~contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
    
    pA=[];
    lbA=[];
    ubA=[];
    
end


if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
    
    pC=log10(InvValues.C);
    lbC=log10(CtrlVar.Cmin)+zeros(size(pC));
    ubC=log10(CtrlVar.Cmax)+zeros(size(pC));
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'c')
    
    pC=InvValues.C;
    lbC=CtrlVar.Cmin+zeros(size(pC));
    ubC=CtrlVar.Cmax+zeros(size(pC));
    
elseif ~contains(lower(CtrlVar.Inverse.InvertFor),'c')
    
    pC=[];
    lbC=[];
    ubC=[];
    
end

p=[pA;pC];
plb=[lbA;lbC];
pub=[ubA;ubC];






end

