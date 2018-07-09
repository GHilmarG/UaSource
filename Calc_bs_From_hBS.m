function [b,s,h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow)

nargoutchk(4,4)
narginchk(7,7)

if ~ isstruct(CtrlVar)
   error('Calc_bs_From_hBSL:InputError','Incorrect inputs.')
end

%% Calculates b, s, and h, consistent with the floating condition.
% [b,s,h,GF]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar)
% calculates b and s from h, S and B and the densities rho and rhow
%
% Note: h is only modified if on input h is smaller than CtrlVar.ThickMin,
%       and CtrlVar.ResetThicknessToMinThickness is true.
%
%
% Step 1: where grounded, b is set equal to B
%         where afloat,   b is calculated from h, rhow and rhow, using the
%         floating condition.
% Step 2: s=b+h
%



% Does the thickness need to be modified?

if CtrlVar.ResetThicknessToMinThickness
    
    h(h<CtrlVar.ThickMin)=CtrlVar.ThickMin;
    %fprintf(CtrlVar.fidlog,' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(indh0),CtrlVar.ThickMin,min(h));
    fprintf(CtrlVar.fidlog,' Setting h(h<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;
end


% Step 1:  
hf=rhow*(S-B)./rho ;

GF.node = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1 if grounded, 0 if afloat

%GF.ele=Nodes2EleMean(MUA.connectivity,GF.node);

bfloat=S-rho.*h/rhow;

b=GF.node.*B + (1-GF.node) .* bfloat ;


% because the grounding line is `smeared out' a bit for a finite CtrlVar.kH
% one can have situations where b<B. For CtrlVar.kH>0.1 this is not really much of an issue
I=b<B ;

if any(I)
    
    if CtrlVar.Report_if_b_less_than_B
        fprintf(CtrlVar.fidlog,' Calc_bs_From_hBS: Found %-i cases where b<B. Setting b>=B.  \n ',numel(find(I))) ;
        
        if CtrlVar.doplots==1
            figure ; plot(MUA.coordinates(I,1)/CtrlVar.PlotXYscale,MUA.coordinates(I,2)/CtrlVar.PlotXYscale,'.') ; axis equal ; title('locations where b<B')
            figure ; plot3(MUA.coordinates(I,1)/CtrlVar.PlotXYscale,MUA.coordinates(I,2)/CtrlVar.PlotXYscale,b(I)-B(I),'.')  ; title('b-B (where negative)')
            
        end
    end
    
    b(I)=B(I); % make sure that lower ice surface is never below bedrock
    
end


% Step 2:
s=b+h;


end