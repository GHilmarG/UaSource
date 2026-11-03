function Code=SlidingLawCode(SlidingLaw)

%% Integer code of a sliding law, from its name or one of its aliases
%
%   Code=SlidingLawCode(CtrlVar.SlidingLaw)
%
%   Code   Name         Aliases
%    1     "Weertman"   "-W-"
%    2     "Budd"       "-B-", "W-N0"
%    3     "Coulomb"    "-C-"
%    4     "Tsai"       "minCW-N0"
%    5     "Cornford"   "rpCW-N0"
%    6     "Umbi"       "rCW-N0"
%    7     "Joughin"    "rCW-v0"
%
% This is the one place where the names and the aliases of the sliding laws are defined (used by BasalDrag.m and BasalDragdC.m).
% The former short aliases "W", "B" and "C" are no longer accepted; use "-W-", "-B-" and "-C-" instead.
%
% Written so that it can be used by MATLAB Coder (character comparisons only, no string arrays).
%
% See also: BasalDrag, BasalDragdC, SlidingLawBeta
%
% (26 Oct 2026)
%%


L=char(SlidingLaw);

if strcmp(L,'Weertman') || strcmp(L,'-W-')
    Code=1;
elseif strcmp(L,'Budd') || strcmp(L,'-B-') || strcmp(L,'W-N0')
    Code=2;
elseif strcmp(L,'Coulomb') || strcmp(L,'-C-')
    Code=3;
elseif strcmp(L,'Tsai') || strcmp(L,'minCW-N0')
    Code=4;
elseif strcmp(L,'Cornford') || strcmp(L,'rpCW-N0')
    Code=5;
elseif strcmp(L,'Umbi') || strcmp(L,'rCW-N0')
    Code=6;
elseif strcmp(L,'Joughin') || strcmp(L,'rCW-v0')
    Code=7;
elseif strcmp(L,'W') || strcmp(L,'B') || strcmp(L,'C')
    error('SlidingLawCode:FormerAlias', ...
        'The sliding-law alias "%s" is no longer used. Use "-W-" (Weertman), "-B-" (Budd) or "-C-" (Coulomb), or the full name.',L)
else
    error('SlidingLawCode:UnknownSlidingLaw','Unknown sliding law "%s".',L)
end

end
