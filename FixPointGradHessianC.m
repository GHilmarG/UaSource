function [dIdC,ddIdCC,UserVar]=FixPointGradHessianC(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo)


%  Fix point estimate of dIdC and ddIdCC


uErr=sqrt(spdiags(Meas.usCov)); vErr=sqrt(spdiags(Meas.vsCov));
usres=(F.ub-Meas.us)./uErr;  vsres=(F.vb-Meas.vs)./vErr;

dIdC=(usres.*F.ub./uErr+vsres.*F.vb./vErr)./F.C;
dIdC=dIdC.*F.GF.node;

ddIdCC=((F.ub./uErr).^2+(F.vb./vErr).^2  )./(F.C.^2);

if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
    
    % ddIddpFP=log(10)^2*F.C.*dIdpFP+(log(10).*F.C).^2 .* ddIddpFP;
    % ddIddpFP=(log(10).*F.C).^2 .* ddIddpFP;
    % dIdpFP=log(10)*F.C.*dIdpFP;
    
    dIdC=log(10)*(usres.*F.ub./uErr+vsres.*F.vb./vErr);
    ddIdCC=(log(10))^2*((F.ub./uErr).^2+(F.vb./vErr).^2 + (CtrlVar.SpeedZero./uErr).^2 );
    
end


dIdC=dIdC.*F.GF.node;

% This should be a FE matrix assembly,   
% ddIdCC=sparse(1:MUA.Nnodes,1:MUA.Nnodes,ddIdCC);

%% Assembly

dIdC=(MUA.M*dIdC)/MUA.Area ; 

[UserVar,ddIdCC]=FE_outer_product(UserVar,CtrlVar,MUA,ddIdCC) ; 
ddIdCC=ddIdCC/MUA.Area; 

ddIdCC=(ddIdCC+ddIdCC')/2 ; 


end