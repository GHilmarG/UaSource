function [UserVar,RunInfo,CtrlVar,MUA]=...
    Remeshing(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,GF,...
    xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened)


if contains(CtrlVar.MeshRefinementMethod,'local','IgnoreCase',true)
    
    MUA=LocalMeshRefinement(CtrlVar,MUA,ElementsToBeRefined,ElementsToBeCoarsened);
    
elseif   contains(CtrlVar.MeshRefinementMethod,'global','IgnoreCase',true)


    [UserVar,RunInfo,MUA]=GlobalRemeshing(UserVar,RunInfo,CtrlVar,MUA,xNod,yNod,EleSizeDesired);
    
else
    fprintf('Incorrect value for CtrlVar.MeshRefinementMethod (%s) \n',CtrlVar.MeshRefinementMethod)
    error('RemeshingBasedOnExplicitErrorEstimate:CaseNotFound','Case not found.')
end



end

