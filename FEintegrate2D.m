function Int=FEintegrate2D(CtrlVar,MUA,f)
%%
%  Int=FEintegrate2D(CtrlVar,MUA,f) calculates the integral of a nodal variable
%  f over each of the elements of the FE mesh.
%
%  CtrlVar can be entered as an empty variable.
%
% Examples: 
% 
% Calculate volume of ice:
%
%   Int=FEintegrate2D([],MUA,h); 
%   TotalIceVolume=sum(Int)
% 
% Calculate grounded area:
%
%   Int=FEintegrate2D([],MUA,GF.node); GroundedArea=sum(Int);
% 
%% 

ndim=2; 

fnod=reshape(f(MUA.connectivity,1),MUA.Nele,MUA.nod);

[points,weights]=sample('triangle',MUA.nip,ndim);

Int=zeros(MUA.Nele,1);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; 
    detJ=MUA.DetJ(:,Iint);
    detJw=detJ*weights(Iint);
    fint=fnod*fun;
    Int=Int+fint.*detJw;

end


end