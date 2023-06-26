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

l=sqrt(2*MUA.EleAreas); 
                       

dtcritical=min(l./(speed+eps));


if all(speed==0)
    dtcritical=nan;
end

end

