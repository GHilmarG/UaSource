function  [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

%  Fully-Implicit solution for uvh (SSTREAM) and h (SSHEET)

narginchk(8,8)
nargoutchk(4,5)

switch lower(CtrlVar.FlowApproximation)
    
    
    case 'sstream'
        
        [UserVar,RunInfo,F1,l1,BCs1]=SSTREAM_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

        
        F1.ud=zeros(MUA.Nnodes,1)  ; F1.vd=zeros(MUA.Nnodes,1); 
        
    case 'ssheet'
        
        
        [UserVar,F1.ud,F1.vd,F1.h,F1.s,F1.GF,l1.h,RunInfo]=SSHEET_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs1,CtrlVar.dt,F1.h,F0.h,F1.S,F1.B,F0.as,F0.ab,F1.as,F1.ab,l1.h,F1.AGlen,F1.n,F1.rho,F1.rhow,F1.g);
     
        
        F1.ub=zeros(MUA.Nnodes,1) ; 
        F1.vb=zeros(MUA.Nnodes,1) ; 
        
    case 'hybrid'
        
        error('uvh2D:CaseNotImplemented','Implicit transient runs not yet implemented for the hybrid flow approximation')
        
    otherwise

        error('what case')
        
end


end


