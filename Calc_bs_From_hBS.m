function [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,coordinates)
    
    
    %% Calculates b, s, and h, consistent with the floating condition.
    % [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar)
    % sets b and b given h, S and B and the densities rho and rhow
    %
    % Note: h is only modified if: 
    %       CtrlVar.ResetThicknessToMinThickness is true,
    %       and, on input h is smaller than CtrlVar.ThickMin
    %
    % Step 1: where grounded, b is set equal to B
    %         where afloat,   b is calculated from h using the floating condition
    % Step 2: s=b+h
    %
    
      
    %if nargin~=7  ; error(' number of input arguments incorrect ' ) ; end
    %if nargout~=3 ; error(' number of output arguments incorrect ' ) ; end
    
    
    if CtrlVar.ResetThicknessToMinThickness==1
        
        h(h<CtrlVar.ThickMin)=CtrlVar.ThickMin;
        %fprintf(CtrlVar.fidlog,' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(indh0),CtrlVar.ThickMin,min(h));
        fprintf(CtrlVar.fidlog,' Setting h(h<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;
    end
    
    %     %%
    %     df=S-rho.*h/rhow;
    %     gf = HeavisideApprox(CtrlVar.kH,B-df,CtrlVar.Hh0);  % 1 if grounded, 0 if afloat
    %
    %     figure ; plot3(coordinates(:,1)/CtrlVar.PlotXYscale,coordinates(:,2)/CtrlVar.PlotXYscale,gf,'.')  ; title('gf')

    
    hf=rhow*(S-B)./rho ;
    gf = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1 if grounded, 0 if afloat
    
    bfloat=S-rho.*h/rhow;
    
    b=gf.*B + (1-gf) .* bfloat ;
    
    
    % because the grounding line is `smeared out' a bit for a finite CtrlVar.kH
    % one can have situations where b<B. For CtrlVar.kH>0.1 this is not really much of an issue
    I=b<B ;
    
    if any(I)
        
        if CtrlVar.Report_if_b_less_than_B
            fprintf(CtrlVar.fidlog,' Calc_bs_From_hBS: Found %-i cases where b<B. Setting b>=B.  \n ',numel(find(I))) ;
            
            if CtrlVar.doplots==1
                figure ; plot(coordinates(I,1)/CtrlVar.PlotXYscale,coordinates(I,2)/CtrlVar.PlotXYscale,'.') ; axis equal ; title('locations where b<B')
                figure ; plot3(coordinates(I,1)/CtrlVar.PlotXYscale,coordinates(I,2)/CtrlVar.PlotXYscale,b(I)-B(I),'.')  ; title('b-B (where negative)')
                
            end
        end
        
        b(I)=B(I); % make sure that lower ice surface is never below bedrock
        
    end
    
    s=b+h;
    
    
end