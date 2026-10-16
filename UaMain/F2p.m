function [p,plb,pub]=F2p(CtrlVar,MUA,F,Meas)


narginchk(4,4)

% p is the vector of the control variables, currently p=[A,B,C]
% with A, B or C here only being nonempty when inverted for,
% 

pA=[];
lbA=[];
ubA=[];

pC=[];
lbC=[];
ubC=[];

pB=[];
lbB=[];
ubB=[];


[isA,isB,isC] = isABC(CtrlVar);

if isA 
    
    pA=log10(F.AGlen);
    
    lbA=log10(F.AGlenmin)+zeros(size(pA));
    ubA=log10(F.AGlenmax)+zeros(size(pA));
        
end


if isC
    
    pC=log10(F.C);
    lbC=log10(F.Cmin)+zeros(size(pC));
    ubC=log10(F.Cmax)+zeros(size(pC));
    
end

if isB

    pB=F.B;

 
    %% Set upper and lower limits for B
    %
    % We have two constraints:
    %
    % 1) B (the bedrock) must be below F.s (the upper glacier surface). Otherwise the ice thickness (F.h) over grounded areas becomes
    % negative. The ice thickness is 
    % 
    %   F.h=F.s-F.b
    % 
    % and where the ice is grounded we have
    %
    % F.B=F.b
    %
    % This condition can be therefore expressed as
    % 
    % B < F.s
    %
    % 2) Additionally, we want to enforce that the ice that was grounded at the beginning of the iteration, never becomes
    % un-grounded. 
    %
    % This can be expressed as 
    %
    %   B>Bstar
    %
    % where
    %
    %   Bstar
    %
    % is the bedrock elevation at flotation.  
    %
    % We only apply B>Bstar where the ice is grounded, i.e. where the nodal grounded/flotation mask, F.GF.node, is greater than 1/2
    %
    % So the condition is
    %
    % B>Bstar where F.GF.node>0.5
    %
    % Similarly, we don't want ice which was afloat to become grounded, i.e. 
    % 
    %
    % B<Bstar where F.GF.node<0.5
    %
    % This second situation is not going to happen if we do not update the bed, B, where the ice is already afloat. 
    %
    % 
    % BAboveFloatationMinimum = 10/CtrlVar.kH ;   % ~10 grounding-line widths
    % 
    % Bstar=(F.s-F.S.*F.rhow./F.rho)./(1-F.rhow./F.rho)+ BAboveFloatationMinimum;  % we must have B > Bstar
    % 
    % GF=F.GF.node>0.5;
    % lbB=nan(MUA.Nnodes,1);
    % lbB(GF)=Bstar(GF) ;      % where grounded, set lower bound just above flotation as based on s, S and densities 
    % BfloatationOffsetBound=100;
    % lbB(~GF)=F.B(~GF)-BfloatationOffsetBound ;  % where afloat, set lower bound to some offset value below the current B. 
    %                                             % I'm not expecting this to matter much as the B where the ice is afloat will not change much during the inversion.
    %                                             % However, because of regularization on B, I might still have some changes in B across the grounding line.

   F.s=Meas.s;
   [F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow);

    % Since I'm now allowing flotation within the domain, I don't have my old lower bounds on B which were there to ensure nothing
    % would go afloat.  But I still need to populate the lbB vector. I guess it might be OK to set it to -inf, but for the time
    % being I simply set to to a save lower limit 
   lbB=-5000 + zeros(MUA.Nnodes,1);

   ubB=F.s-CtrlVar.ThickMin ; % This is Bmax, we must have B < Bmax = s - thickmin

   if ~isempty(lbB)
       ubB=max(lbB,ubB) ; % make sure ubB >= lbB
   end

   % UaPlots(CtrlVar,MUA,F,lbB,FigureTitle="lbB")
   % UaPlots(CtrlVar,MUA,F,ubB,FigureTitle="ubB")
    % UaPlots(CtrlVar,MUA,F,ubB-lbB,FigureTitle="ubB-lbB") ; CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
    % UaPlots(CtrlVar,MUA,F,"-B-",FigureTitle="B")

end


p=[pA;pB;pC];
plb=[lbA;lbB;lbC];
pub=[ubA;ubB;ubC];

% make sure is feasible
p=kk_proj(p,pub,plb) ;


%% Box transformation
%
% Eliminate the box constraints by the logistic change of variables, p -> u, see BoxTransform.m . After this the problem
% is genuinely unconstrained, so setting plb and pub empty is now a correct description rather than a temporary measure.
%
% This has to come BEFORE the Cholesky mapping below. The Cholesky mapping mixes components, so a box in p would become
% a general polytope in the mapped variable, which is why box constraints could not be used with it. Applied after the
% box transformation there is no box left to mix up, and the two compose cleanly.
%
% The bounds are taken from MUA rather than from the plb and pub just computed, so that exactly the same box is used
% here, in p2F and in JGH, and for the whole run. See the note in InvertForModelParameters.m .

if CtrlVar.Inverse.BoxTransform

    if ~isfield(MUA,'BoxTransform') || isempty(MUA.BoxTransform)
        error('F2p:NoStoredBox',...
            ['CtrlVar.Inverse.BoxTransform is true but MUA.BoxTransform has not been set.\n',...
            'It must be built once, before the inversion starts, in InvertForModelParameters.m .'])
    end

    p=BoxTransform("forward",p,MUA.BoxTransform) ;

    plb=[];
    pub=[];

end


if CtrlVar.Inverse.CholeskyMappingOfCostFunctionAndGradient

    % If the box transformation was applied above, plb and pub are already empty and the mapping acts on a genuinely
    % unconstrained variable, which is the case it was designed for. Otherwise the old behaviour applies: the box
    % constraints cannot be carried through the mapping and are discarded, with a note.
    if ~isempty(plb) || ~isempty(pub)
        plb=[];
        pub=[];
        fprintf("Note: For Cholesky mapping of cost function and gradient, box constraints can not be used.\n")
        fprintf("      Box constraints are now eliminated. \n")
        fprintf("      No box constraints on any of the inverted fields are used.\n")
        fprintf("      Set CtrlVar.Inverse.BoxTransform=true to eliminate them by a change of variables instead,\n")
        fprintf("      which retains them rather than discarding them.\n")
    end

    % p=MUA.RG*p;
    p=MUA.RG*(MUA.PRG'*p);

end


end

