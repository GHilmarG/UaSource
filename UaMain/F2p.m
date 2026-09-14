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


[isA,isB,isC] = isABC(CtrlVar);

if isA 
    
    pA=log10(F.AGlen);
    
    lbA=log10(F.AGlenmin)+zeros(size(pA));
    ubA=log10(F.AGlenmax)+zeros(size(pA));
        
end


if isC
    
    pC=log10(F.C);
    lbC=log10(F.Cmin)+zeros(size(pC));
    ubC=log10(F.Cmax)+zeros(size(pC));
    
end

if isB

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


if CtrlVar.Inverse.CholeskyMappingOfCostFunctionAndGradient


    % p=MUA.RG*p;
    p=MUA.RG*(MUA.PRG'*p);

    plb=[];
    pub=[];
    fprintf("Note: For Cholesky mapping of cost function and gradient, box constraints can not be used.\n")
    fprintf("      Box constraints are now eliminated. \n")
    fprintf("      No box constraints on any of the inverted fields are used.\n")


end


end

