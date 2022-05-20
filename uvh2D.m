function  [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

%  Fully-Implicit solution for uvh (SSTREAM) and h (SSHEET)

narginchk(8,8)
nargoutchk(4,5)

switch lower(CtrlVar.FlowApproximation)
    
    
    case {"sstream","sstream-rho"}
        
        [UserVar,RunInfo,F1,l1,BCs1]=SSTREAM_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

        
        F1.ud=zeros(MUA.Nnodes,1)  ; F1.vd=zeros(MUA.Nnodes,1); 
        
    case {'ssheet','hybrid'}

        [UserVar,RunInfo,F1,l1,BCs1]=SSHEET_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

    case "uvhprescribed"


        [UserVar,RunInfo,F1,l1,BCs1]=uvhPrescibed(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

    otherwise

        error('what case')
        
end


end


