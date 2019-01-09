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
        
        
    case 'B'
        
        % tested and works
        % express geometrical variables in terms of p
        
        F.B=p ; 
        F.h= F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p)  ;   % h = s - b

        %         bfloat=F.S - F.rho.*F.h /F.rhow;
%         
%         F.b=F.GF.node.*p + (1-F.GF.node) .* bfloat ;
%            =F.GF.node.*p + (1-F.GF.node) .* (F.S - F.rho.*F.h /F.rhow) 
%            =F.GF.node.*p + (1-F.GF.node) .* (F.S - F.rho.*(F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p))  ;  
          F.b=F.GF.node.*p + (1-F.GF.node) .* (F.S - F.rho.*( F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p)   )./F.rhow);
          
          % b = GF  p + (1-GF) (S-rho (h0 (1-GF) + GF (s0-p) /rhow)
          % b = GF  p + (1-GF) (S-rho (h0 (1-GF) + GF (s0-p) /rhow)
          % db/dp = GF + (1-GF) rho GF/rhow
        
        %         dB/dp = 1
        %         dh/dp = -GF.node
        %         db/dp = GF.node + (1-GF.node) (0 - rho dhdp/rhow)
        %               = GF.node + (1-GF.node) (0 -- rho GF.node/rhow)
        %               = GF.node + (1-GF.node) (+ rho GF.node/rhow)
        %               = GF.node + rho (1-GF.node) GF.node/ rhow
        %
        
        [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow);
        
         % F.h= F.hInit.*(1-F.GF.node)  + F.GF.node.* (F.sInit - p)  ;  % because GF has changed
        
    case 'b'
        
        F.h=F.s-p ;
        F.B=p.*F.GF.node+(1-F.GF.node).*F.BInit;
        
        
        %  bfloat=F.S - F.rho.*(F.s-p) /F.rhow;
        %  dbfloat/dp= F.rho./F.rhow
        %
        % F.h=F.GF.node.*(F.s-F.B)+(1-F.GF.node).*F.hInit ;
        
        [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow);
        
        
        
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