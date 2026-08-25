






function [KHess_qp]=Hess_qp(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC)

narginchk(13,13)



 [KFCu,KFCv]=FCuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;

 [KFAu,KFAv]=FAuv(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y) ;


KFpq=[KFAu KFAv ; ...
      KFCu KFCv ] ;

xi=[KdudA KdudC ;...
    KdvdA KdvdC] ; 

K=KFpq*xi; 

KHess_qp=K+K' ; 



end