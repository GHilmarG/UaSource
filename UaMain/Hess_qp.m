


function [KHess_qp]=Hess_qp(CtrlVar,MUA,F,Psi_x,Psi_y,KdudA,KdvdA,KdhdA,KdudB,KdvdB,KdhdB,KdudC,KdvdC,KdhdC)



[Psid2FdAdu,Psid2FdAdv]=Psi_d2FdAdq_xi(CtrlVar,MUA,F,Psi_x,Psi_y);
[Psid2FdCdu,Psid2FdCdv]=Psi_d2FdCdq_xi(CtrlVar,MUA,F,Psi_x,Psi_y);



%
% Fpqξ=⟨Ψ, δpq2​F[ϕl​,ξ,m​]⟩
% ⟨Ψ, δqp2​F[ξ,l​,ϕm​]⟩=(Fpqξ)T
%

Psid2Fdpq=[Psid2FdAdu Psid2FdAdv ; ...
           Psid2FdCdu Psid2FdCdv ] ;

xi=[KdudA KdudC ...
    KdvdA KdvdC] ; 

K=Psid2Fdpq*xi; 

KHess_qp=K+K' ; 



end