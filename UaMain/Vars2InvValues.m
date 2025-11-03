
function InvValues=Vars2InvValues(CtrlVar,F,InvValues,J,dJdp,JGHouts,RunInfo,dJdpTest)

NA=numel(F.AGlen);
NB=numel(F.B);
NC=numel(F.C);


% always return as inverse final values, the corresponding F fields (as suggested by Camilla)
InvValues.AGlen=F.AGlen;
InvValues.B=F.B;
InvValues.C=F.C;

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
            
            
            InvValues.dJdBTest=dJdpTest(1:NB);
            InvValues.dJdCTest=dJdpTest(NB+1:end);
            
        case 'ABC'
            
            
            InvValues.dJdAGlenTest=dJdpTest(1:NA);
            InvValues.dJdBTest=dJdpTest(NA+1:NA+NB);
            InvValues.dJdCTest=dJdpTest(NA+NB+1:end);
            
            
            
        otherwise
            
            error('case error')
    end
    
 
end


end

