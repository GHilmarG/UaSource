






function [KHess_qp]=Hess_qp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC)

narginchk(16,16)

% [KFAu,KFAv]=Psi_d2FdAdq_xi(CtrlVar,MUA,F,Psi_x,Psi_y);
% [KFCu,KFCv]=Psi_d2FdCdq_xi(CtrlVar,MUA,F,Psi_x,Psi_y);

 [KFCu,KFCv]=FCuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;

 [KFAu,KFAv]=FAuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;


KFpq=[KFAu KFAv ; ...
      KFCu KFCv ] ;




xi=[KdudA KdudC ;...
    KdvdA KdvdC] ; 


K=KFpq*xi; 

KHess_qp=K+K' ; 



end