function Int=FEintergrate2D(CtrlVar,MUA,f)

%%
%  Int=FEintergrate2D(CtrlVar,MUA,f)
%  calculates the integral of f over each of the elements of the FE mesh given in the structure MUA
%  CtrlVar can be entered as an empty variable
%
%
%

ndim=2; 

fnod=reshape(f(MUA.connectivity,1),MUA.Nele,MUA.nod);

[points,weights]=sample('triangle',MUA.nip,ndim);

Int=zeros(MUA.Nele,1);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ')
        detJ=MUA.DetJ(:,Iint);
    else
        [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    
    detJw=detJ*weights(Iint);
    fint=fnod*fun;
    
    Int=Int+fint.*detJw;

end


end