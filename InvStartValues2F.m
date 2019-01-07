function F=InvStartValues2F(CtrlVar,F,InvStartValues,Priors)

if  contains(CtrlVar.Inverse.InvertForField,'A')
    
    if isempty(InvStartValues.AGlen)
        
        fprintf('InvStartValues.AGlen can not be left empty when inverting for A.\n')
        error('InvStartValues2F:InvStartValues.AGlen')
        
    end
    
    if isempty(Priors.AGlenmin)
        
        fprintf('Priors.AGlenmin can not be left empty when inverting for A.\n')
        error('InvStartValues2F:Priors.AGlenminEmpty')
        
    end
    
    
    F.AGlen=InvStartValues.AGlen ;
    F.AGlenmin=Priors.AGlenmin;
    F.AGlenmax=Priors.AGlenmax;
    
end

if  contains(CtrlVar.Inverse.InvertForField,'b')
    
    if isempty(InvStartValues.b)
        
        fprintf('InvStartValues.b can not be left empty when inverting for b or B.\n')
        error('InvStartValues2F:InvStartValues.b')
        
    end
    
    F.b=InvStartValues.b ;
    F.bmin=Priors.bmin;
    F.bmax=Priors.bmax;
    
end

if  contains(CtrlVar.Inverse.InvertForField,'B')
    
    if isempty(InvStartValues.B)
        
        fprintf('InvStartValues.B can not be left empty when inverting for b or B.\n')
        error('InvStartValues2F:InvStartValues.b')
        
    end
    
    F.B=InvStartValues.B ;
    F.Bmin=Priors.Bmin;
    F.Bmax=Priors.Bmax;
    
end



if  contains(CtrlVar.Inverse.InvertForField,'C')
    
    if isempty(InvStartValues.C)
        
        fprintf('InvStartValues.C can not be left empty when inverting for C.\n')
        error('InvStartValues2F:InvStartValues.C')
        
    end
    
    F.C=InvStartValues.C ;
    F.Cmin=Priors.Cmin;
    F.Cmax=Priors.Cmax;
    
end

end

