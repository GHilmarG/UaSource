function Mask=CalcMeshMask(CtrlVar,MUA,NodalValue,Level)

%%
%
%%

narginchk(3,6)
nargoutchk(1,4)

if isempty(NodalValue)
    return
end


EleValue=Nodes2EleMean(MUA.connectivity,NodalValue);

CtrlVar.GLthreshold=Level ; % need this in GLgeometry
Mask.Level=CtrlVar.GLthreshold; 
[Mask.Geo,Mask.NodesOn,Mask.ElementsOn]=GLgeometry(MUA.connectivity,MUA.coordinates,NodalValue,CtrlVar);

Mask.ElementsIn=   (EleValue>Mask.Level)  & ~Mask.ElementsOn;
Mask.ElementsOut= (EleValue<Mask.Level)  & ~Mask.ElementsOn;


Mask.NodesIn=false(MUA.Nnodes,1);
Mask.NodesOut=false(MUA.Nnodes,1);

Mask.NodesIn(MUA.connectivity(Mask.ElementsIn,:))=true;
Mask.NodesOut(MUA.connectivity(Mask.ElementsOut,:))=true;


Mask.NodesOut=Mask.NodesOut & ~Mask.NodesOn;
Mask.NodesIn=Mask.NodesIn & ~Mask.NodesOn;


return



end