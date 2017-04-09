function [UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=GetDesiredEleSize(UserVar,CtrlVar,MUA,F,GF,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators)


narginchk(11,11)
nargoutchk(4,4)

N=nargout('DefineDesiredEleSize');

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
        
end

end

