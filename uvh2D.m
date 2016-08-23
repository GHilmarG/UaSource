function [UserVar,ub1,vb1,ud1,vd1,h1,lambdauv1,lambdah1,RunInfo]=...
    uvh2D(UserVar,CtrlVar,MUA,BCs,dt,h0,S,B,ub0,vb0,ud0,vd0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,dudt,dvdt,lambdauv,lambdah,...
    AGlen,C,n,m,alpha,rho,rhow,g)

%  Fully-Implicit solution for uvh (SSTREAM) and h (SSHEET)


nargoutchk(9,9)

switch lower(CtrlVar.FlowApproximation)
    
    
    case 'sstream'
        
        [UserVar,ub1,vb1,h1,lambdauv1,lambdah1,RunInfo]=...
            SSTREAM_TransientImplicit(UserVar,CtrlVar,MUA,BCs,dt,h0,S,B,ub0,vb0,ub1,vb1,h1,as0,ab0,as1,ab1,dudt,dvdt,lambdauv,lambdah,...
            AGlen,C,n,m,alpha,rho,rhow,g);
        ud1=ud1*0 ; vd1=vd1*0 ; 
        
    case 'ssheet'
        
        
        [UserVar,ud1,vd1,h1,s1,lambdah1,RunInfo]=SSHEET_TransientImplicit(UserVar,CtrlVar,MUA,BCs,dt,h1,h0,S,B,as0,ab0,as1,ab1,lambdah,AGlen,n,rho,rhow,g);
        lambdauv1=lambdauv;
        
        ub1=ub1*0 ; vb1=vb1*0;
        
    case 'hybrid'
        
        error('uvh2D:CaseNotImplemented','Implicit transient runs not yet implemented for the hybrid flow approximation')
        
    otherwise

        error('what case')
        
end


end


