


function [g,J] = JuvLSQ(x,RunInfo,CtrlVar,MUA,F,BCs)


% g : gradient
% J : Jacobian

% RunInfo=UaRunInfo();
% l=[]; 
% [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv([],RunInfo,CtrlVar,MUA,BCs,F,l);
%
F.ub=x(1:MUA.Nnodes);
F.vb=x(MUA.Nnodes+1:end);

CtrlVar.uvAssembly.ZeroFields=false ; 
if nargout ==1
    CtrlVar.uvMatrixAssembly.Ronly=true;

    [RunInfo,Ruv]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F);
    Kuv=[]; 



else

    CtrlVar.uvMatrixAssembly.Ronly=false;
    [RunInfo,Ruv,Kuv]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F,BCs); 


end

g=Ruv;
J=Kuv; 

end















