function d=Draft(S,B,b,CtrlVar)
    
    d = HeavisideApprox((1/CtrlVar.GroundingLineTransition.lmin),S-B,CtrlVar.Hh0).*(S-b);  % draft
    
end
