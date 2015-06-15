function [h1,lambdah]=SSS2dPrognostic(dt,h0,u0,v0,du0dt,dv0dt,a0,da0dt,u1,v1,a1,da1dt,du1dt,dv1dt,coordinates,connectivity,Boundary,nip,L,Lrhs,lambdah,Itime,CtrlVar)

% semi-implicit forward integration.
% implicit with respect to thickness (h)
% explicit with respect to velocity

VectorVersion=1; nonVectorVersion=0;


switch lower(CtrlVar.FlowApproximation)
    
    
    case 'sstream';
        
        
        if nonVectorVersion==1
            [h1,lambdah]=...
                Nexh2DSparse(dt,h0,u0,v0,a0,u1,v1,a1,coordinates,connectivity,Nnodes,Nele,nip,nod,L,Lrhs,lambdah,CtlrVar.theta,Itime,CtrlVar);
        end
        
        if VectorVersion==1
            
            if CtrlVar.TG3==1
                [h1,lambdah]=NexthTG3in2D(dt,h0,u0,v0,du0dt,dv0dt,a0,da0dt,u1,v1,a1,da1dt,du1dt,dv1dt,coordinates,connectivity,Boundary,nip,L,Lrhs,lambdah,CtrlVar);
            else
                [h1,lambdah]=Nexh2DSparseVector(dt,h0,u0,v0,a0,u1,v1,a1,coordinates,connectivity,nip,L,Lrhs,lambdah,Itime,CtrlVar);
            end
        end
        
    otherwise
        
        error('Semi-implicit time integration only implemented for SSTREAM and not SSHEET or Hybrid formulations \n')
        
end

end
