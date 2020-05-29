function dtcritical=CalcCFLdt2D(UserVar,RunInfo,CtrlVar,MUA,F)

%
%  CFL   : dt < l/v
%
%


switch lower(CtrlVar.FlowApproximation)
    
    
    case 'sstream'
        
        u=F.ub;
        v=F.vb;
        
        
    case 'ssheet'
        
        u=F.ud ;
        v=F.vd ;
        
        
    case 'hybrid'
        
        u=F.ub+F.ud ;
        v=F.vb+F.vd ;
        
        
    otherwise
        
        error('what case')
        
end

speed=sqrt(u.*u+v.*v) ;
speed=Nodes2EleMean(MUA.connectivity,speed); 

l=sqrt(MUA.EleAreas);  % there is a factor of two here that I'm not including, this will reduced dt 
                       % and provide a margin of safety.

dtcritical=min(l./(speed+eps));

end

