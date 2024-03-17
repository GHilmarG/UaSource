function [b,s,h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow)

nargoutchk(2,4)
narginchk(7,7)

if ~ isstruct(CtrlVar)
    error('Calc_bs_From_hBSL:InputError','Incorrect inputs.')
end

%% Calculates b, s, and h, consistent with the floating condition.
%
%    [b,s,h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow)
%
%
% Note: h is only modified if on input h is smaller than CtrlVar.ThickMin,
%       and CtrlVar.ResetThicknessToMinThickness is true.
%
%   MUA is only needed for plotting purposes, can be left empty.
%
% Step 1: where grounded, b is set equal to B
%         where afloat,   b is calculated from h, rhow and rhow, using the
%         floating condition.
% Step 2: s=b+h
%
%%


% Does the thickness need to be modified?

if CtrlVar.ResetThicknessToMinThickness

    indh0=h<CtrlVar.ThickMin;

    if any(indh0)
        h(indh0)=CtrlVar.ThickMin;

        if CtrlVar.InfoLevelThickMin>=1
            fprintf(' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(find(indh0)),CtrlVar.ThickMin,min(h));
            fprintf(' Setting h(h<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;

            if CtrlVar.InfoLevelThickMin>=10 && CtrlVar.doplots
                fig=FindOrCreateFigure('ThickMin');
                PlotMuaMesh(CtrlVar,MUA)
                hold on
                plot(MUA.coordinates(indh0,1)/CtrlVar.PlotXYscale,MUA.coordinates(indh0,2)/CtrlVar.PlotXYscale,'or')
            end
        end
    end
end


% Step 1:
%if nargin<=8
hf=rhow.*(S-B)./rho ;
GF.node = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1 if grounded, 0 if afloat
%end

%GF.ele=Nodes2EleMean(MUA.connectivity,GF.node);

bfloat=S-rho.*h./rhow;

b=GF.node.*B + (1-GF.node) .* bfloat ;


if CtrlVar.Enforce_bAboveB
    % because the grounding line is `smeared out' a bit for a finite CtrlVar.kH one can
    % have situations where b<B. For CtrlVar.kH>0.1 this is not really much of an issue
    % and anyhow this is a direct consequence of using a smooth step function (explained in detail in UaCompendium)
    I=b<B ;

    if any(I)

        if ( CtrlVar.InfoLevelThickMin>=1)  && ~isempty(MUA)
            fprintf(CtrlVar.fidlog,' Calc_bs_From_hBS: Found %-i cases where b<B. Setting b>=B.  \n ',numel(find(I))) ;

            if CtrlVar.doplots==1
                fig=FindOrCreateFigure('b<B');
                PlotMuaMesh(CtrlVar,MUA) ;
                hold on ;
                plot(MUA.coordinates(I,1)/CtrlVar.PlotXYscale,MUA.coordinates(I,2)/CtrlVar.PlotXYscale,'.r') ; title('locations where b<B')
                PlotGroundingLines(CtrlVar,MUA,GF);

            end
        end

        b(I)=B(I); % make sure that lower ice surface is never below bedrock

    end
end

% Step 2:
s=b+h;

if CtrlVar.GroundingFloatingMaskContains=="GF nodes and strickly afloat/grounded nodes and elements"
    GF=IceSheetIceShelves(CtrlVar,MUA,GF);
end



end