
function InvValues=Vars2InvValues(CtrlVar,F,InvValues,J,dJdp,JGHouts,RunInfo,dJdpTest)

NA=numel(F.AGlen);
NB=numel(F.B);
NC=numel(F.C);


% always return as inverse final values, the corresponding F fields (as suggested by Camilla)
InvValues.AGlen=F.AGlen;
InvValues.B=F.B;
InvValues.C=F.C;

InvValues.J=J;
InvValues.I=JGHouts.I;
InvValues.R=JGHouts.R;
InvValues.RAGlen=[];
InvValues.RC=[];

InvValues.RCa=[];
InvValues.RCs=[];
InvValues.RAa=[];
InvValues.RAs=[];

if isfield(JGHouts,'Psi_x') && isfield(JGHouts,'Psi_y')
    InvValues.uAdjoint=JGHouts.Psi_x;
    InvValues.vAdjoint=JGHouts.Psi_y;
else
    InvValues.uAdjoint=[];
    InvValues.vAdjoint=[];
end


%% Gradients

InvValues.dJdp=dJdp;

InvValues.dIdp=JGHouts.dIdp;
InvValues.dRdp=JGHouts.dRdp;
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

