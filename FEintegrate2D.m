





function Int=FEintegrate2D(CtrlVar,MUA,f)
%%
%  Int=FEintegrate2D(CtrlVar,MUA,f) calculates the integral of a nodal or integration-point variable
%  f over each of the elements of the FE mesh.
%
%  CtrlVar can be entered as an empty variable.
%
% *Examples:*
% 
% Ice volume:
%
%   Int=FEintegrate2D([],MUA,h); 
%   TotalIceVolume=sum(Int)
% 
% Grounded area:
%
%   Int=FEintegrate2D([],MUA,GF.node); GroundedArea=sum(Int);
% 
% Integrated surface mass balance over grounded areas (i.e. upstream of grounding lines):
%
%  M=sum(FEintegrate2D(CtrlVar,MUA,F.as.*F.GF.node)) ;
%
%
%% 

ndim=2; 



%%
%
% 
% 
% [points,weights]=sample('triangle',MUA.nip,ndim);

Int=zeros(MUA.Nele,1);


if isempty(MUA.DetJ)
    MUA=UpdateMUA(CtrlVar,MUA) ;
end


if numel(f)==MUA.Nnodes  % nodal variable

    fnod=reshape(f(MUA.connectivity,1),MUA.Nele,MUA.nod);
    for Iint=1:MUA.nip

        fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
        detJ=MUA.DetJ(:,Iint);
        detJw=detJ*MUA.weights(Iint);
        fint=fnod*fun;
        Int=Int+fint.*detJw;

    end

elseif size(f,1)==MUA.Nele && size(f,2) == MUA.nip  % integration point variable

    for Iint=1:MUA.nip

        
        detJ=MUA.DetJ(:,Iint);
        detJw=detJ*MUA.weights(Iint);
        fint=f(:,Iint);
        Int=Int+fint.*detJw;

    end

end