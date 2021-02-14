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
    
    % F.b=InvStartValues.b ;
    fprintf(' Note:  In a B inversion, the upper ice surface (s) is set equal to Meas.s .\n')
    fprintf('        The user must ensure that Meas.s is consistent with floatation. \n')
    fprintf('        The lower surface (b) and the grounding-line will initially be calculated using Meas.s.\n') 
    F.s=Meas.s ;
    
    % First calculate b from s, B and S given rho and rhow.
    [F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow,F.GF); 
    % Now again calculate b, s from h, S and B to ensure full consistency with 
    % the rest of the code.  This might results in s changing and being different
    % from Meas.s
    [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow) ;
    sDiff=norm(F.s-Meas.s);
    fprintf(' After having calculated b from s, and then again s from h=s-h using floation we find norm(F.s-Meas.s)=%g \n',sDiff)
    fprintf(' Now we set Meas.s=F.s using the new F.s calculated from floation.\n')
    Meas.s=F.s ; % 
end



if  contains(CtrlVar.Inverse.InvertForField,'C')
    
    if isempty(InvStartValues.C)
        
        fprintf('InvStartValues.C can not be left empty when inverting for C.\n')
        error('InvStartValues2F:InvStartValues.C')
        
    end
    
    F.C=InvStartValues.C ;
    F.n=InvStartValues.n;
    F.muk=InvStartValues.muk ;
    F.q=InvStartValues.q;
    
    
    F.Cmin=Priors.Cmin;
    F.Cmax=Priors.Cmax;
    
end

end

