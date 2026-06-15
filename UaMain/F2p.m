function [p,plb,pub]=F2p(CtrlVar,MUA,F)


narginchk(3,3)

% p is the vector of the control variables, currently p=[A,B,C]
% with A, B or C here only being nonempty when inverted for,
% 

pA=[];
lbA=[];
ubA=[];

pC=[];
lbC=[];
ubC=[];

pB=[];
lbB=[];
ubB=[];

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

    pB=F.B;

    % lbB=F.Bmin+zeros(size(pB));
    % ubB=F.Bmax+zeros(size(pB));



    % Where initially grounded, make sure ice never goes afloat
    BAboveFloatationMinimum=10;
    Bstar=(F.s-F.S.*F.rhow./F.rho)./(1-F.rhow./F.rho)+ BAboveFloatationMinimum;  % this is Bmin, we must have B > Bmin

    GF=F.GF.node>0.5;
    lbB=nan(MUA.Nnodes,1);
    lbB(GF)=Bstar(GF) ;      % where grounded, set lower bound just above flotation as based on s, S and densities 
    lbB(~GF)=F.B(~GF)-100 ;  % where afloat, set lower bound to some small value, although this should not really have an impact on retrieved B



    % set


    % ensure that min ice thickness is not violated
    %ubB=[];
    ubB=F.s-CtrlVar.ThickMin ; % This is Bmax, we must have B < Bmax

    ubB=max(lbB,ubB) ; % make sure ubB >= lbB

    % UaPlots(CtrlVar,MUA,F,lbB,FigureTitle="lbB")
    % UaPlots(CtrlVar,MUA,F,ubB,FigureTitle="ubB")
    % UaPlots(CtrlVar,MUA,F,ubB-lbB,FigureTitle="ubB-lbB") ; CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
    % UaPlots(CtrlVar,MUA,F,"-B-",FigureTitle="B")

end


p=[pA;pB;pC];
plb=[lbA;lbB;lbC];
pub=[ubA;ubB;ubC];

% make sure is feasible
p=kk_proj(p,pub,plb) ;

end

