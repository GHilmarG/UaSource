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
        
        
        %F.b=F.GFInit.node.*p+(1-F.GFInit.node).*F.bInit ;
        %F.B=F.GFInit.node.*p+(1-F.GFInit.node).*F.BInit ;
        
        %F.b=F.GF.node.*p+(1-F.GF.node).*F.bInit ;
        %F.B=F.GF.node.*p+(1-F.GF.node).*F.BInit ;
        
        bOld=F.b;    
        BOld=F.B;

        if contains(lower(CtrlVar.Inverse.InvertFor),'-B-')
            
            F.b=F.GF.node.*p+(1-F.GF.node).*F.b ;
            
        else
            
            F.b=p;
            
        end
        
        Ver=4 ;
        
        switch Ver
            
            case 1
                
                % [----------- This produces almost correct results:
                dh=-(F.b-bOld); % h smaller if b raised
                F.h=F.h+dh;
                %F.B=min([F.b F.B],[],2);  % not a good idea
                F.B=F.GF.node.*F.b+(1-F.GF.node).*BOld ; 
                F.s=F.b+F.h;
                % ------------]
            case 2

                
                F.B=F.GF.node.*F.b+(1-F.GF.node).*BOld ;
      
                db=(F.b-bOld); % h smaller if b raised
                
                dh=-F.GF.node.*db-(1-F.GF.node).*db.*F.rhow./F.rho ;
                
                nk=0;
                for k=1:nk
                    
                    hf=F.rhow*(F.S-F.B)./F.rho ;
                    F.GF.node = HeavisideApprox(CtrlVar.kH,F.h+dh-hf,CtrlVar.Hh0);
                    dh=-F.GF.node.*db-(1-F.GF.node).*db.*F.rhow./F.rho ;
                    F.B=F.GF.node.*F.b+(1-F.GF.node).*BOld ;
                    
                end
                
                F.h=F.h+dh;
                F.s=F.b+F.h;
                
                
            case 3
                % [----------- This produces almost correct results too
                
                F.h=F.s-F.b;
                F.B=F.GF.node.*F.b+(1-F.GF.node).*BOld ;
                
            case 4
                
                 F.h=F.s-F.b;
                 F.B=F.GF.node.*F.b+(1-F.GF.node).*BOld ; % this is needed 
%                 hf=F.rhow*(F.S-F.B)./F.rho ;
%                 F.GF.node = HeavisideApprox(CtrlVar.kH,F.h-hf,CtrlVar.Hh0);
%                 F.B=F.GF.node.*F.b+(1-F.GF.node).*BOld ;
        end
        
        
        % In the assembly I do not use F.b! only s and h 
        % I recalculate at He and deltas
        % dbdx=dsdx-dhdx
        
        %[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow);
        
        
        
        % Need to change F.s across the grounding line!
        % but do it in such a way that the F.b does not change!
        % Its OK if h changes again a bit, as I'm modifying h to be consitent
        % with changes in B
        
%         hf=F.rhow*(F.S-F.B)./F.rho ;
%         F.GF.node = HeavisideApprox(CtrlVar.kH,F.h-hf,CtrlVar.Hh0);
%         F.s=F.GF.node.*(F.h+F.B)+ (1-F.GF.node).*(F.S+(1-F.rho/F.rhow).*F.h); 
%         F.h=F.s-F.b;
%         
        %[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,[],F.h,F.S,F.B,F.rho,F.rhow);
        
     
        
        
        
        
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