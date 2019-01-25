function F=InvStartValues2F(CtrlVar,MUA,F,InvStartValues,Priors,Meas)

narginchk(6,6)

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
    error('sfdaf')
    
end

if    contains(CtrlVar.Inverse.InvertForField,'B')
    
    
    if isempty(InvStartValues.B)
        
        fprintf('InvStartValues.B can not be left empty when inverting for b or B.\n')
        error('InvStartValues2F:InvStartValues.b')
        
    end
    
    F.B=InvStartValues.B ;
    F.Bmin=Priors.Bmin;
    F.Bmax=Priors.Bmax;
    
    
    
    F.s=Meas.s ; % note that since I'm not inverting for s, I must keep s fixed, this of course
    % may not be possible over the floating areas, consider calculating F.b over the floating areas
    % from F.s using the floating relationship.
    
    F.b=Calc_b_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow,F.GF);
    
    % However, this now creates the possibility that b<B in places downstream of the groundign line.
    % What to do?
    % Either: Modify the floating mask by keeping h fixed and recalculating s, b, B, and GF.
    %         First modify B so that B<b everywhere, then recalcuate s, b, B and GF
    
    F.B=min(F.b,F.B); 
    [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.s-F.b,F.S,F.B,F.rho,F.rhow);
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

