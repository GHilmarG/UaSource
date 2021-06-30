function Int=FEintegrate2Dint(CtrlVar,MUA,fint)

%%
%  Int=FEintegrate2Dint(CtrlVar,MUA,fint) 
%%

Int=zeros(MUA.Nele,1);

for Iint=1:MUA.nip

    detJ=MUA.DetJ(:,Iint);
    detJw=detJ*MUA.weights(Iint);
    Int=Int+fint.*detJw;

end


end