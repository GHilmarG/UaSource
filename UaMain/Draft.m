function d=Draft(S,B,b,CtrlVar)
    
    d = HeavisideApprox(CtrlVar.kH,S-B,CtrlVar.Hh0).*(S-b);  % draft
    
end
