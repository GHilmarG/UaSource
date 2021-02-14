function [dIdA,ddIdAA,UserVar]=FixPointGradHessianA(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo)


%  Fix point estimate of dIdC and ddIdCC

dIdA=[];

uErr=sqrt(spdiags(Meas.usCov)); vErr=sqrt(spdiags(Meas.vsCov));

ddIdAA=((F.ub./uErr).^2+(F.vb./vErr).^2  )./(F.AGlen.^2);



if contains(lower(CtrlVar.Inverse.InvertFor),'loga')
    
    ddIdAA=(log(10))^2*((F.ub./uErr).^2+(F.vb./vErr).^2 + (CtrlVar.SpeedZero./uErr).^2 );
end


% This should be a FE matrix assembly,   
% ddIdAC=sparse(1:MUA.Nnodes,1:MUA.Nnodes,ddIdAC);

%% Assembly



[UserVar,ddIdAA]=FE_outer_product(UserVar,CtrlVar,MUA,ddIdAA) ; 
ddIdAA=ddIdAA/MUA.Area; 
ddIdAA=(ddIdAA+ddIdAA')/2 ; 


end