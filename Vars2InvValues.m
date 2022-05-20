
function InvValues=Vars2InvValues(CtrlVar,F,InvValues,J,dJdp,JGHouts,RunInfo,dJdpTest)

NA=numel(F.AGlen);
Nb=numel(F.b);
NC=numel(F.C);


if contains(CtrlVar.Inverse.InvertForField,'A')
    
    InvValues.AGlen=F.AGlen;
    
end

if contains(CtrlVar.Inverse.InvertForField,'B')
    

    InvValues.B=F.B;
    
    % This should be consistent with:
    %[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
end




if contains(CtrlVar.Inverse.InvertForField,'C')
    
    InvValues.C=F.C;
    
end



InvValues.J=J;
InvValues.I=JGHouts.MisfitOuts.I;
InvValues.R=JGHouts.RegOuts.R;
InvValues.RAGlen=JGHouts.RegOuts.RAGlen;
InvValues.RC=JGHouts.RegOuts.RC;

InvValues.RCa=JGHouts.RegOuts.RCa;
InvValues.RCs=JGHouts.RegOuts.RCs;
InvValues.RAa=JGHouts.RegOuts.RAa;
InvValues.RAs=JGHouts.RegOuts.RAs;

if isprop(InvValues,'uAdjoint')
    InvValues.uAdjoint=JGHouts.MisfitOuts.uAdjoint;
    InvValues.vAdjoint=JGHouts.MisfitOuts.vAdjoint;
end


%% Gradients

InvValues.dJdp=dJdp;

InvValues.dIdp=JGHouts.dIdp;
InvValues.dRdp=JGHouts.dRdp;

InvValues.dJdAGlen=JGHouts.MisfitOuts.dIdAGlen+JGHouts.RegOuts.dRdAGlen;
InvValues.dJdC=JGHouts.MisfitOuts.dIdC+JGHouts.RegOuts.dRdC;
InvValues.dJdB=NaN; % 
InvValues.dJdB=JGHouts.MisfitOuts.dIdB+JGHouts.RegOuts.dRdB;

%% These are of less interest, but can be added
%InvFinalValues.dIdAGlen=JGHouts.MisfitOuts.dIdAGlen;
%InvFinalValues.dIdC=JGHouts.MisfitOuts.dIdC;

%InvFinalValues.dRdAGlen=JGHouts.RegOuts.dRdAGlen;
%InvFinalValues.dRdC=JGHouts.RegOuts.dRdC;


InvValues.SearchStepSize=RunInfo.Inverse.StepSize(end);

%%
if ~isempty(dJdpTest)
    
    InvValues.dJdAGlenTest=[];
    InvValues.dJdBTest=[];
    InvValues.dJdCTest=[];
    
    InvValues.dJdpTest=dJdpTest;
    
    switch  CtrlVar.Inverse.InvertForField
        
        case 'A'
            
            InvValues.dJdAGlenTest=dJdpTest;

            
        case 'B'
            
            InvValues.dJdBTest=dJdpTest;
            
        case 'C'
            
            InvValues.dJdCTest=dJdpTest;
  
            
        case 'AC'
            
            InvValues.dJdAGlenTest=dJdpTest(1:NA);
            InvValues.dJdCTest=dJdpTest(NA+1:end);
            
        case 'BC'
            
            
            InvValues.dJdBTest=dJdpTest(1:Nb);
            InvValues.dJdCTest=dJdpTest(Nb+1:end);
            
        case 'ABC'
            
            
            InvValues.dJdAGlenTest=dJdpTest(1:NA);
            InvValues.dJdBTest=dJdpTest(NA+1:NA+Nb);
            InvValues.dJdCTest=dJdpTest(NA+Nb+1:end);
            
            
            
        otherwise
            
            error('case error')
    end
    
 
end


end

