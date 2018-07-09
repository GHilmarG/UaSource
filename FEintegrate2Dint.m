function Int=FEintegrate2Dint(CtrlVar,MUA,fint)

%%
%  Int=FEintegrate2Dint(CtrlVar,MUA,fint) 
%%

ndim=2; 

%fnod=reshape(f(MUA.connectivity,1),MUA.Nele,MUA.nod);

[points,weights]=sample('triangle',MUA.nip,ndim);

Int=zeros(MUA.Nele,1);

for Iint=1:MUA.nip
    
    %fun=shape_fun(Iint,ndim,MUA.nod,points) ; 
    detJ=MUA.DetJ(:,Iint);
    detJw=detJ*weights(Iint);
    %fint=fnod*fun;
    Int=Int+fint.*detJw;

end


end