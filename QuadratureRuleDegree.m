function Degree=QuadratureRuleDegree(CtrlVar)

switch CtrlVar.TriNodes
    
    case 3 
        
        Degree=4;
        
        
    case 6
        
        Degree=8;
        
        
    case 10
        
        Degree=8 ;
        
    otherwise
        
        error('Ua:CaseNotFound','Which case?')
        
end




end