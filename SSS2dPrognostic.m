function [UserVar,RunInfo,h1,l]=SSS2dPrognostic(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)

%  function [UserVar,RunInfo,h1,l]=SSS2dPrognostic(UserVar,RunInfo,CtrlVar,MUA,BCs,l,h0,ub0,vb0,dub0dt,dvb0dt,a0,da0dt,ub1,vb1,a1,da1dt,dub1dt,dvb1dt)

nargoutchk(3,4)


% semi-implicit forward integration.
% implicit with respect to thickness (h)
% explicit with respect to velocity



MLC=BCs2MLC(CtrlVar,MUA,BCs);
Lh=MLC.hL ; Lhrhs=MLC.hRhs ;

if numel(l.h)==numel(Lhrhs)
    lambdah=l.h;   % TO DO/to do: consider checking that l.h has indeed the right dimensions. Here lambdah is only an initial guess for lambdah when solved using an iterative solver.
else
    lambdah=Lhrhs*0;
end

dt=CtrlVar.dt;

switch lower(CtrlVar.FlowApproximation)
    
    
    case 'sstream'
        
        switch CtrlVar.uvhSemiImplicitTimeSteppingMethod
            
            % for CtrlVar.TG3=0 both do the same, but Next2DSparseVector does not calculate
            % the TG3 terms, where as NextTG3in2D always does, even if they are not needed.
            case "TG3"
                % This always includes the TG3 terms
                [h1,lambdah]=NexthTG3in2D(dt,h0,ub0,vb0,dub0dt,dvb0dt,a0,da0dt,ub1,vb1,a1,da1dt,dub1dt,dvb1dt,MUA.coordinates,MUA.connectivity,MUA.Boundary,MUA.nip,Lh,Lhrhs,lambdah,CtrlVar);
                l.h=lambdah;
            case "Galerkin"
                % This is based on NexthTG3in2D but does not include the TG3 terms
                % This used to be the default approach (until early 2020), but it does not
                % include SUPG terms
                [h1,lambdah]=Nexh2DSparseVector(dt,h0,ub0,vb0,a0,ub1,vb1,a1,MUA.coordinates,MUA.connectivity,MUA.nip,Lh,Lhrhs,lambdah,CtrlVar);
                l.h=lambdah;
            case "SUPG"
                
                % kappa=zeros(MUA.Nnodes,1);
                % [UserVar,h1,lambdah]=TracerConservationEquation(UserVar,CtrlVar,MUA,dt,h0,ub0,vb0,a0,ub1,vb1,a1,kappa,BCs);
                
                [UserVar,RunInfo,h1,l]=MassContinuityEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
                
        end
        
        
    otherwise
        
        error('Semi-implicit time integration only implemented for SSTREAM and not SSHEET or Hybrid formulations \n')
        
end



end
