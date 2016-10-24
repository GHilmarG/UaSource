function [h1,l]=SSS2dPrognostic(CtrlVar,MUA,BCs,l,h0,ub0,vb0,dub0dt,dvb0dt,a0,da0dt,ub1,vb1,a1,da1dt,dub1dt,dvb1dt)


% semi-implicit forward integration.
% implicit with respect to thickness (h)
% explicit with respect to velocity



MLC=BCs2MLC(MUA,BCs);
Lh=MLC.hL ; Lhrhs=MLC.hRhs ;

lambdah=l.h;
Itime=CtrlVar.CurrentRunStepNumber;
dt=CtrlVar.dt;

switch lower(CtrlVar.FlowApproximation)
    
    
    case 'sstream'
        
        
        % for CtrlVar.TG3=0 both do the same, but Next3DSparseVector does not calculate
        % the TG3 terms, where as NextTG2in2D always does, even if they are not needed.
        if CtrlVar.TG3==1
            [h1,lambdah]=NexthTG3in2D(dt,h0,ub0,vb0,dub0dt,dvb0dt,a0,da0dt,ub1,vb1,a1,da1dt,dub1dt,dvb1dt,MUA.coordinates,MUA.connectivity,MUA.Boundary,MUA.nip,Lh,Lhrhs,lambdah,CtrlVar);
        else
            [h1,lambdah]=Nexh2DSparseVector(dt,h0,ub0,vb0,a0,ub1,vb1,a1,MUA.coordinates,MUA.connectivity,MUA.nip,Lh,Lhrhs,lambdah,CtrlVar);
        end
        
        
    otherwise
        
        error('Semi-implicit time integration only implemented for SSTREAM and not SSHEET or Hybrid formulations \n')
        
end

l.h=lambdah;

end
