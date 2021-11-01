function [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=GetDesiredEleSize(UserVar,CtrlVar,MUA,F,GF,xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators)


persistent nCalls

narginchk(11,11)
nargoutchk(4,4)

InputFile="DefineDesiredEleSize.m"; TestIfInputFileInWorkingDirectory(InputFile) ; 


if isempty(nCalls)
    nCalls=0;
end
nCalls=nCalls+1;


N=nargout('DefineDesiredEleSize');

if contains(lower(CtrlVar.MeshRefinementMethod),"local")
    x=MUA.xEle ; y=MUA.yEle;
else
    x=MUA.coordinates(:,1) ;
    y=MUA.coordinates(:,2) ;
end

InputEleSizeDesired=EleSizeDesired;
InputElementsToBeRefined=ElementsToBeRefined;

switch N
    
    case 3
        
        
        fprintf('\n Note:  DefineDesiredEleSize only has three output arguments, but now allows for four.\n')
        fprintf('          The new fourth output argument is ''ElementsToBeCoarsened'' \n\n')
        
        if nargin(InputFile)>5
            
            [UserVar,EleSizeDesired,ElementsToBeRefined]=...
                DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,F.s,F.b,F.S,F.B,F.rho,F.rhow,F.ub,F.vb,F.ud,F.vd,GF,NodalErrorIndicators);
        else
            [UserVar,EleSizeDesired,ElementsToBeRefined]=...
                DefineDesiredEleSize(UserVar,CtrlVar,MUA,F,x,y,EleSizeDesired,ElementsToBeRefined,NodalErrorIndicators);
        end
        
        ElementsToBeCoarsened=[];
        
        
    case 4
        
        if nargin(InputFile)>5
            
            [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
                DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,F.s,F.b,F.S,F.B,F.rho,F.rhow,F.ub,F.vb,F.ud,F.vd,GF,NodalErrorIndicators);
        else
            
             [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
                DefineDesiredEleSize(UserVar,CtrlVar,MUA,F,x,y,EleSizeDesired,ElementsToBeRefined,NodalErrorIndicators);
            
        end
        
    otherwise
        
        
        error('Ua:GetDesiredEleSize','DefineDesiredEleSize.m must return either 3 or 4 output arguments')
        
        
end

if contains(CtrlVar.MeshRefinementMethod,'local','IgnoreCase',true)
    
    if ~isequal(InputEleSizeDesired,EleSizeDesired)
        if nCalls==1
            fprintf(['Note: The mesh refinement methods is %s. \n' ...
                'When using this local mesh-refinement method in combination with ``DefineDesiredEleSize'' \n',...
                'you only need to specify as an output which elements are to be refined and/or coarsended. \n',...
                'You can specify ``EleSizeDesired'' as well, but it will be ignored.\n'],...
                CtrlVar.MeshRefinementMethod)
        end
    end
    
    % just some basic check here, not covering all possibilies, just the most likely input mistakes

    if islogical(ElementsToBeRefined)
        
        if numel(ElementsToBeRefined) ~= MUA.Nele
            
            fprintf("Number of elements in the logical list ElementsToBeRefined (%i) not equal to number of elements in mesh (%i) \n ",...
                numel(ElementsToBeRefined),MUA.Nele)
            error("Ua:IncorectUserInput","Incorrect DefineDesiredEleSize user input")
        end
        
    end

   if islogical(ElementsToBeCoarsened)
        
        if numel(ElementsToBeCoarsened) ~= MUA.Nele
            
            fprintf("Number of elements in the logical list ElementsToBeCoarsened (%i) not equal to number of elements in mesh (%i) \n ",...
                numel(ElementsToBeCoarsened),MUA.Nele)
            error("Ua:IncorectUserInput","Incorrect DefineDesiredEleSize user input")
        end
        
    end




elseif contains(CtrlVar.MeshRefinementMethod,'global','IgnoreCase',true)
    
    if ~isequal(InputElementsToBeRefined,ElementsToBeRefined)
        if nCalls==1
            fprintf(['Note: The mesh refinement methods is %s. \n' ...
                'When using this global mesh-refinement method in combination with ``DefineDesiredEleSize'' \n',...
                'you only need to specify as an output the variable ``EleSizeDesired''.\n',...
                'You can specify which elements are to be refined and/or coarsended but this will be ignored. \n'],...
                CtrlVar.MeshRefinementMethod)
        end
    end
    
    
else
    fprintf('Incorrect value for CtrlVar.MeshRefinementMethod (%s) \n',CtrlVar.MeshRefinementMethod)
    error('RemeshingBasedOnExplicitErrorEstimate:CaseNotFound','Case not found.')
end




end

