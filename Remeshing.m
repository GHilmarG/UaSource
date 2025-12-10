function [UserVar,RunInfo,CtrlVar,MUAnew]=...
    Remeshing(UserVar,RunInfo,CtrlVar,MUAold,BCs,F,l,GF,...
    xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened)


if contains(CtrlVar.MeshRefinementMethod,'local','IgnoreCase',true)
    
    [MUAnew,RunInfo]=LocalMeshRefinement(CtrlVar,RunInfo,MUAold,ElementsToBeRefined,ElementsToBeCoarsened);
    
elseif contains(CtrlVar.MeshRefinementMethod,'global','IgnoreCase',true)
    
    F.x=MUAold.coordinates(:,1); F.y=MUAold.coordinates(:,2);
    [UserVar,RunInfo,MUAnew]=GlobalRemeshing(UserVar,RunInfo,CtrlVar,MUAold,xNod,yNod,EleSizeDesired,F);
    
else
    fprintf('Incorrect value for CtrlVar.MeshRefinementMethod (%s) \n',CtrlVar.MeshRefinementMethod)
    error('RemeshingBasedOnExplicitErrorEstimate:CaseNotFound','Case not found.')
end





end

