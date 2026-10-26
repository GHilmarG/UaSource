function [R,ok,missing]=MakeKernelParamsFlex(CtrlVar)
% Run-time parameters for the flexible compiled uvh kernel: a fixed set of fields with fixed types (double scalars and char rows).
% ok=false if a required field is missing in CtrlVar, or an option is not supported (the vectorised assembly is then used instead).
%
% (25 Oct 2026) The lists of fields are parsed only once (persistent), as this is called in each uvh assembly.

persistent P
if isempty(P)
    Num=["AGlenmin","AGlenminWidth","Cmin","CminWidth","Czero","EpsZero","Hh0","Nzero","SpeedZero","ThickMin","etaZero","kH","theta","dt", ...
     "LevelSetMethodMassBalanceFeedbackCoeffCubic","LevelSetMethodMassBalanceFeedbackCoeffLin", ...
     "ThicknessPenaltyMassBalanceFeedbackCoeffCubic","ThicknessPenaltyMassBalanceFeedbackCoeffLin","ThicknessPenaltyMassBalanceFeedbackCoeffQuad", ...
     "ThicknessPenaltyMassBalanceFeedbackExponential.K","ThicknessPenaltyMassBalanceFeedbackExponential.l", ...
     "ThicknessPenaltyMassBalanceFeedbackSoftPlus.K","ThicknessPenaltyMassBalanceFeedbackSoftPlus.deltaAbs","ThicknessPenaltyMassBalanceFeedbackSoftPlus.l","ThicknessPenaltyMassBalanceFeedbackSoftPlus.lRel", ...
     "SUPG.beta0","uvh.SUPG.tauMultiplier","Hydrology.gammaw","ThicknessBarrierMassBalanceFeedbackCoeffLog", ...
     "Enforce_bAboveB","IncludeMelangeModelPhysics","LevelSetMethod","LevelSetMethodAutomaticallyApplyMassBalanceFeedback", ...
     "OnlyCalcBasalDragAndEffectiveViscosity","ResetThicknessToMinThickness","ThicknessBarrier","ThicknessPenalty","ThicknessPenaltyNodal","Inverse.dFuvdClambda"];
    Chr=["SlidingLaw","ThicknessPenaltyMassBalanceFeedbackFunction","uvh.SUPG.tau"];
    % fields that may be missing in CtrlVar, with the value the code assumes if they are missing
    Opt=struct("ThicknessBarrier",0,"ThicknessBarrierMassBalanceFeedbackCoeffLog",0,"OnlyCalcBasalDragAndEffectiveViscosity",0,"ThicknessPenaltyNodal",0,"Inverse_dFuvdClambda",0,"Hydrology_gammaw",NaN);
    P.Num=cellfun(@(p) cellstr(split(p,".")).',cellstr(Num),'UniformOutput',false) ; P.NumName=Num ;
    P.Chr=cellfun(@(p) cellstr(split(p,".")).',cellstr(Chr),'UniformOutput',false) ; P.ChrName=Chr ;
    P.Opt=Opt ; P.OptKey=strrep(Num,".","_") ;
end

R=struct(); ok=true; missing=strings(0);
for k=1:numel(P.Num)
    [in,v]=getpath(CtrlVar,P.Num{k});
    if ~in
        if isfield(P.Opt,P.OptKey(k)) ; v=P.Opt.(P.OptKey(k)) ; else ; ok=false ; missing(end+1)=P.NumName(k) ; v=NaN ; end %#ok<AGROW>
    end
    if isempty(v) ; v=NaN ; end
    if ~isscalar(v) || ~(isnumeric(v) || islogical(v)) ; ok=false ; missing(end+1)=P.NumName(k) ; v=NaN ; end %#ok<AGROW>
    R=setpath(R,P.Num{k},double(v));
end
for k=1:numel(P.Chr)
    [in,v]=getpath(CtrlVar,P.Chr{k});
    if ~in ; ok=false ; missing(end+1)=P.ChrName(k) ; v='' ; end %#ok<AGROW>
    R=setpath(R,P.Chr{k},char(v));
end
% display-only fields: plotting is not part of the compiled code
R.doplots=0; R.PlotSUPGparameter=0; R.InfoLevelThickMin=0;
% integer codes for the string options, as in the forward branch of BasalDrag.m and in SUPGtau.m (0: not supported, use the vectorised assembly)
Laws={{'Weertman','W'},{'Budd','W-N0'},{'Coulomb','C'},{'Tsai','minCW-N0'},{'rpCW-N0','Cornford'},{'rCW-N0','Umbi'},{'Joughin','rCW-v0'}};
R.SlidingLawCode=0; for k=1:numel(Laws) ; if any(strcmp(R.SlidingLaw,Laws{k})) ; R.SlidingLawCode=k ; end ; end
if R.SlidingLawCode==0 ; ok=false ; missing(end+1)="SlidingLaw="+string(R.SlidingLaw)+" (not supported)" ; end
R.uvh.SUPG.tauCode=find(strcmp(R.uvh.SUPG.tau,{'tau1','tau2','taus','taut','tau0'}));
if isempty(R.uvh.SUPG.tauCode) ; R.uvh.SUPG.tauCode=0 ; ok=false ; missing(end+1)="uvh.SUPG.tau="+string(R.uvh.SUPG.tau)+" (not supported)" ; end
% the inverse-run branch of BasalDrag (dFuvdC) is not part of the compiled kernel
if R.Inverse.dFuvdClambda~=0 ; ok=false ; missing(end+1)="Inverse.dFuvdClambda (inverse branch not compiled)" ; end
end

function [in,v]=getpath(S,p)
in=true; v=S;
for j=1:numel(p)
    if isstruct(v) && isfield(v,p{j}) ; v=v.(p{j}) ; else ; in=false ; v=[] ; return ; end
end
end

function S=setpath(S,p,v)
switch numel(p)
    case 1 ; S.(p{1})=v ;
    case 2 ; S.(p{1}).(p{2})=v ;
    otherwise ; S.(p{1}).(p{2}).(p{3})=v ;
end
end
