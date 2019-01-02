function InvValues=p2InvValues(CtrlVar,p,InvValues,NA,NC)


if contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')  % AC
    
    
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
        InvValues.AGlen=10.^p(1:NA);
    else
        InvValues.AGlen=p(1:NA);
    end
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        InvValues.C=10.^p(NA+1:end);
    else
        InvValues.C=p(NA+1:end);
    end
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'aglen')   % A
    
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
        InvValues.AGlen=10.^p;
    else
        InvValues.AGlen=p;
    end
    
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'c')  % C
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        InvValues.C=10.^p;
    else
        InvValues.C=p;
    end
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'-b-')  % C
    
    
    %  ad gera
    InvValues.b=p;
    % I=GF.node>0.5; %only change b and B where grounded
    % InvValues.b(I)=p(I); % this does change the thickness
    % InvValues.B(I)=InvValues.b(I); % now change B where grounded
    % F.h=F.s-F.b;

%    [F.b,F.s,F.h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
    
else
    fprintf(' CtrlVar.Inverse.InvertFor=%s \n',CtrlVar.Inverse.InvertFor)
    fprintf(' CtrlVar.Inverse.InvertFor does not have expected value.\n')
    error('InverForModelParameters:incorrect inputs')
end


end