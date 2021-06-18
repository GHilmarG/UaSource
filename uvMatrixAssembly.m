function [Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(CtrlVar,MUA,F)

%
% Ruv=Tint-Fext;
% Tint   : internal nodal forces
% Fint   : external nodal forces

narginchk(3,5)
nargoutchk(1,5)

switch lower(CtrlVar.FlowApproximation)
    
    case "sstream"
        
        [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F) ;
        
    case "sstream-rhotest"
        
        
        [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAMrhoTest(CtrlVar,MUA,F) ;
        
    otherwise
        
        error("UA:uvMatrixAssemblyCaseNotFound","Case Not Found")
        
end


end



