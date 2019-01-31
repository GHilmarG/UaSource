function [p,plb,pub]=F2p(CtrlVar,MUA,F)


narginchk(3,3)

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
    
    pA=log10(F.AGlen);
    
    lbA=log10(F.AGlenmin)+zeros(size(pA));
    ubA=log10(F.AGlenmax)+zeros(size(pA));
    
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'-aglen-')
    
    pA=F.AGlen;
    lbA=F.AGlenmin+zeros(size(pA));
    ubA=F.AGlenmax+zeros(size(pA));
    
end


if contains(lower(CtrlVar.Inverse.InvertFor),'-logc-')
    
    pC=log10(F.C);
    lbC=log10(F.Cmin)+zeros(size(pC));
    ubC=log10(F.Cmax)+zeros(size(pC));
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'-c-')
    
    pC=F.C;
    lbC=F.Cmin+zeros(size(pC));
    ubC=F.Cmax+zeros(size(pC));
    
end

if contains(CtrlVar.Inverse.InvertFor,'-b-')
    
    error('fdsa')
 
    
end


if contains(CtrlVar.Inverse.InvertFor,'-B-')
    
    pb=F.B;
    lbb=F.Bmin+zeros(size(pb));
    ubb=F.Bmax+zeros(size(pb));
    
end


if CtrlVar.Inverse.pPreMultiplier=="M"
    Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);
    if ~isempty(pA)
        pA=(MUA.M*pA)/Area;
        lbA=(MUA.M*lbA)/Area;
        ubA=(MUA.M*ubA)/Area;
    end
    if ~isempty(pb)
        pb=(MUA.M*pb)/Area;
        lbb=(MUA.M*lbb)/Area;
        ubb=(MUA.M*ubb)/Area;
    end
    if ~isempty(pC)
        pC=(MUA.M*pC)/Area;
        lbC=(MUA.M*lbC)/Area;
        ubC=(MUA.M*ubC)/Area;
    end
end

p=[pA;pb;pC];
plb=[lbA;lbb;lbC];
pub=[ubA;ubb;ubC];


end

