function F=p2F(CtrlVar,p,F)

NA=numel(F.AGlen);
Nb=numel(F.b);
NC=numel(F.C);

switch  CtrlVar.Inverse.InvertForField
    
    case 'A'
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            F.AGlen=10.^p(1:NA);
        else
            F.AGlen=p(1:NA);
        end
        
    case 'b'
   
         
        I=F.GF.node>0.5; %only change b and B where grounded
        F.b(I)=p(I);     % this does change the thickness
        F.B(I)=F.b(I);   % now change B where grounded
        F.h=F.s-F.b;

        % [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow); % should not be needed
        
    case 'C'
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            F.C=10.^p;
        else
            F.C=p;
        end
        
    case 'Ab'
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            F.AGlen=10.^p(1:NA);
        else
            F.AGlen=p(1:NA);
        end
        
              
        I=find(F.GF.node>0.5); %only change b and B where grounded
        F.b(I)=p(I+NA);     % this does change the thickness
        F.B(I)=F.b(I);   % now change B where grounded
        F.h=F.s-F.b;
        
       % F.b=p(NA+1:end);
        
        
    case 'AC'
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            F.AGlen=10.^p(1:NA);
        else
            F.AGlen=p(1:NA);
        end
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            F.C=10.^p(NA+1:end);
        else
            F.C=p(NA+1:end);
        end
        
    case 'bC'
        
                   
        I=find(F.GF.node>0.5); %only change b and B where grounded
        F.b(I)=p(I);     % this does change the thickness
        F.B(I)=F.b(I);   % now change B where grounded
        F.h=F.s-F.b;
        
        %F.b=p(1:Nb);
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            F.C=10.^p(Nb+1:end);
        else
            F.C=p(Nb+1:end);
        end
        
        
    case 'AbC'
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            F.AGlen=10.^p(1:NA);
        else
            F.AGlen=p(1:NA);
        end
        
        I=find(F.GF.node>0.5); %only change b and B where grounded
        F.b(I)=p(NA+I);     % this does change the thickness
        F.B(I)=F.b(I);   % now change B where grounded
        F.h=F.s-F.b;
        
        
        %F.b=p(NA+1:NA+Nb);
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            F.C=10.^p(NA+Nb+1:end);
        else
            F.C=p(NA+Nb+1:end);
        end
        
    otherwise
        
        error('p2F: which case?')
end


end