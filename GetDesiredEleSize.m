function [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=GetDesiredEleSize(UserVar,CtrlVar,MUA,F,GF,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators)


narginchk(11,11)
nargoutchk(4,4)

N=nargout('DefineDesiredEleSize');


InputEleSizeDesired=EleSizeDesired;
InputElementsToBeRefined=ElementsToBeRefined;

switch N
    
    case 3
        
        
        fprintf('\n Note:  DefineDesiredEleSize only has three output arguments, but now allows for four.\n')
        fprintf('          The new fourth output argument is ''ElementsToBeCoarsened'' \n\n')
        
        
        [UserVar,EleSizeDesired,ElementsToBeRefined]=...
            DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,F.s,F.b,F.S,F.B,F.rho,F.rhow,F.ub,F.vb,F.ud,F.vd,GF,NodalErrorIndicators);
        
        ElementsToBeCoarsened=[];
        
        
    case 4
        
        [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
            DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,F.s,F.b,F.S,F.B,F.rho,F.rhow,F.ub,F.vb,F.ud,F.vd,GF,NodalErrorIndicators);
        
    otherwise
        
        
        error('Ua:GetDesiredEleSize','DefineDesiredEleSize.m must return either 3 or 4 output arguments')
        
        
end

if contains(CtrlVar.MeshRefinementMethod,'local','IgnoreCase',true)
    
    if ~isequal(InputEleSizeDesired,EleSizeDesired)
        fprintf(['Note: The mesh refinement methods is %s. \n' ...
            'When using this local mesh-refinement method in combination with ``DefineDesiredEleSize'' \n',...
            'you only need to specify as an output which elements are to be refined and/or coarsended. \n',...
            'You can specify ``EleSizeDesired'' as well, but it will be ignored.\n'],...
            CtrlVar.MeshRefinementMethod)
    end
    
    
elseif contains(CtrlVar.MeshRefinementMethod,'global','IgnoreCase',true)
    
    if ~isequal(InputElementsToBeRefined,ElementsToBeRefined)
        fprintf(['Note: The mesh refinement methods is %s. \n' ...
            'When using this global mesh-refinement method in combination with ``DefineDesiredEleSize'' \n',...
            'you only need to specify as an output the variable ``EleSizeDesired''.\n',...
            'You can specify which elements are to be refined and/or coarsended but this will be ignored. \n'],...
            CtrlVar.MeshRefinementMethod)
    end

    
else
    fprintf('Incorrect value for CtrlVar.MeshRefinementMethod (%s) \n',CtrlVar.MeshRefinementMethod)
    error('RemeshingBasedOnExplicitErrorEstimate:CaseNotFound','Case not found.')
end




end

