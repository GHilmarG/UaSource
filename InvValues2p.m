function [p,plb,pub]=InvValues2p(CtrlVar,InvValues)

% p is the vector of the control variables, currenty p=[A,b,C]
% with A, b or C here only being nonempty when inverted for,
% This mapping between A, b and C into the control variable is done by
% InveValues2o

pA=[];
lbA=[];
ubA=[];

pC=[];
lbC=[];
ubC=[];

pb=[];
lbb=[];
ubb=[];

if contains(lower(CtrlVar.Inverse.InvertFor),'-logaglen-')
    
    pA=log10(InvValues.AGlen);
    lbA=log10(CtrlVar.AGlenmin)+zeros(size(pA));
    ubA=log10(CtrlVar.AGlenmax)+zeros(size(pA));
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'-aglen-')
    
    pA=InvValues.AGlen;
    lbA=CtrlVar.AGlenmin+zeros(size(pA));
    ubA=CtrlVar.AGlenmax+zeros(size(pA));
    
end


if contains(lower(CtrlVar.Inverse.InvertFor),'-logc-')
    
    pC=log10(InvValues.C);
    lbC=log10(CtrlVar.Cmin)+zeros(size(pC));
    ubC=log10(CtrlVar.Cmax)+zeros(size(pC));
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'-c-')
    
    pC=InvValues.C;
    lbC=CtrlVar.Cmin+zeros(size(pC));
    ubC=CtrlVar.Cmax+zeros(size(pC));
    
end



if contains(lower(CtrlVar.Inverse.InvertFor),'-b-')
    
    pb=InvValues.b;
    lbb=-1e10+zeros(size(pb));
    ubb=1e10+zeros(size(pb));
end


p=[pA;pb;pC];
plb=[lbA;lbb;lbC];
pub=[ubA;ubb;ubC];


end

