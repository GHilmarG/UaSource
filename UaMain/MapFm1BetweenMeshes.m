function [RunInfo,Fm1new]=MapFm1BetweenMeshes(CtrlVar,RunInfo,MUAold,MUAnew,Fm1old)
%%
% [RunInfo,Fm1new]=MapFm1BetweenMeshes(CtrlVar,RunInfo,MUAold,MUAnew,Fm1old)
%
% Maps the rates of the previous time step (Fm1, see UpdateFtimeDerivatives.m) from the mesh MUAold onto the mesh MUAnew. (9 Oct 2026)
%
% Only the rates dhdt, dubdt, dvbdt, duddt and dvddt, and the time step over which they were calculated (dtRates), are needed in Fm1.
% The rates are interpolated in the same way as the rates of F in MapFbetweenMeshes.m. Nodes of the new mesh outside of the old mesh,
% and rates that are not available (empty or of the wrong size), are set to NaN. The explicit estimate then falls back, node by node,
% to linear extrapolation, see ExplicitEstimationUsingBackwardDifferences.m.
%
%%

Fm1new=UaFields;
Fm1new.dtRates=Fm1old.dtRates;

Names=["dhdt","dubdt","dvbdt","duddt","dvddt"];
isAvailable=false(numel(Names),1);
for k=1:numel(Names)
    isAvailable(k)=numel(Fm1old.(Names(k)))==MUAold.Nnodes;
end

for k=find(~isAvailable)'
    Fm1new.(Names(k))=NaN(MUAnew.Nnodes,1);
end

if any(isAvailable)
    In=cell(1,sum(isAvailable));  Out=cell(1,sum(isAvailable));
    iA=find(isAvailable);
    for k=1:numel(iA)
        In{k}=Fm1old.(Names(iA(k)));
    end
    [RunInfo,Out{:}]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,NaN(1,numel(iA)),In{:});
    for k=1:numel(iA)
        Fm1new.(Names(iA(k)))=Out{k};
    end
end

end
