function nClipped=ReportClipping(CtrlVar,Name,x,lower,upper,Options)

%% Reports values of a field that are outside its lower and upper limits
%
%   nClipped=ReportClipping(CtrlVar,Name,x,lower,upper)
%   nClipped=ReportClipping(CtrlVar,Name,x,lower,upper,LowerName="Cmin",UpperName="Cmax",Inclusive=false,Log=true,Note="",Where="")
%
% Called before the values are set to the limits (eg with kk_proj), so that the user is told that this happens. If
% CtrlVar.InfoLevel>=1 (in an inverse run: CtrlVar.InfoLevel>=1 or CtrlVar.InfoLevelInverse>=1, as inverse runs are often done with
% InfoLevel=0), prints the number and fraction of the values below the lower and above the upper limit, and how far beyond
% the limits they are: in orders of magnitude (Log=true, the default, for positive fields such as C and AGlen), or as a difference
% (Log=false, eg for B, in m, and whenever values or limits are not positive).
%
% A report is printed once for a given caller, field and number of values beyond the limits: in a transient run the slipperiness
% and A are defined at every time step, and the same report would otherwise be repeated at every step. With CtrlVar.InfoLevel>=10
% every report is printed. The record of printed reports is cleared at the start of each run (ClearPersistentUaVariables.m).
%
% Inputs:
%
%   Name          : name of the field, eg "C", "AGlen", "B"
%   x             : the values, before they are set to the limits
%   lower, upper  : the limits, scalars or of the same size as x, or empty (no limit)
%
% Options:
%
%   LowerName, UpperName : names of the limits used in the report (default: Name+"min" and Name+"max")
%   Inclusive     : also count values equal to the limits (default false). For values that may already have been set to the
%                   limits, eg the start values of an inversion.
%   Log           : distance beyond the limits in orders of magnitude (default true) or as a difference (false)
%   Note          : text printed after the report, if any values are beyond the limits
%   Where         : name at the start of the report (default: the calling function)
%
% Output:
%
%   nClipped      : [number of values below the lower limit, number above the upper limit]
%
% Example:
%
%   ReportClipping(CtrlVar,"C",C,CtrlVar.Cmin,CtrlVar.Cmax) ;
%   C=kk_proj(C,CtrlVar.Cmax,CtrlVar.Cmin) ;
%
% (20 Nov 2026)
%%

arguments
    CtrlVar struct
    Name (1,1) string
    x {mustBeNumeric}
    lower {mustBeNumeric}
    upper {mustBeNumeric}
    Options.LowerName (1,1) string = Name+"min"
    Options.UpperName (1,1) string = Name+"max"
    Options.Inclusive (1,1) logical = false
    Options.Log (1,1) logical = true
    Options.Note (1,1) string = ""
    Options.Where (1,1) string = ""
end

persistent Reported

nClipped=[0 0] ;
if isempty(x) ; return ; end

x=x(:) ; N=numel(x) ;
lower=lower(:) ; upper=upper(:) ;

if isempty(lower)
    iL=false(N,1) ;
elseif Options.Inclusive
    iL=x<=lower ;
else
    iL=x<lower ;
end

if isempty(upper)
    iU=false(N,1) ;
elseif Options.Inclusive
    iU=x>=upper ;
else
    iU=x>upper ;
end

nL=nnz(iL) ; nU=nnz(iU) ; nClipped=[nL nU] ;

Level=0 ;
if isfield(CtrlVar,"InfoLevel") ; Level=CtrlVar.InfoLevel ; end
if isfield(CtrlVar,"InverseRun") && CtrlVar.InverseRun && isfield(CtrlVar,"InfoLevelInverse") ; Level=max(Level,CtrlVar.InfoLevelInverse) ; end

if nL+nU==0 || Level<1
    return
end

Where=Options.Where ;
if Where==""
    st=dbstack(1) ;
    if isempty(st) ; Where="ReportClipping" ; else ; Where=string(st(1).name) ; end
end

Key=Where+":"+Name+":"+nL+":"+nU ;
if isempty(Reported) ; Reported=strings(0,1) ; end
if any(Reported==Key) && Level<10
    return
end
Reported(end+1,1)=Key ;

if Options.Inclusive ; Below="are at or below" ; Above="are at or above" ; else ; Below="are below" ; Above="are above" ; end

if nL>0
    Text=LimitReport(x,lower,iL,-1,Options.LowerName,Options.Log) ;
    fprintf("%s: %i of %i %s values (%.3g%%) %s %s and are set to %s.%s\n",Where,nL,N,Name,100*nL/N,Below,LimitText(lower,Options.LowerName),Options.LowerName,Text) ;
end

if nU>0
    Text=LimitReport(x,upper,iU,+1,Options.UpperName,Options.Log) ;
    fprintf("%s: %i of %i %s values (%.3g%%) %s %s and are set to %s.%s\n",Where,nU,N,Name,100*nU/N,Above,LimitText(upper,Options.UpperName),Options.UpperName,Text) ;
end

if Options.Note~=""
    fprintf("%s\n",Options.Note) ;
end

end

%%
function T=LimitText(Limit,LimitName)

if isscalar(Limit)
    T=sprintf("%s=%g",LimitName,Limit) ;
else
    T=sprintf("their limit %s",LimitName) ;
end

end

%%
function T=LimitReport(x,Limit,I,Sign,LimitName,Log)

% how far beyond the limit the values are: in orders of magnitude (if Log and values and limits are positive), otherwise as a
% difference

xI=x(I) ;
if isscalar(Limit) ; LI=Limit+0*xI ; else ; LI=Limit(I) ; end

if all(xI==LI)
    T=" (All of these values are on the limit.)" ;
elseif Log && all(xI>0) && all(LI>0)
    d=Sign*(log10(xI)-log10(LI)) ;
    [dmax,k]=max(d) ;
    T=sprintf(" The value furthest beyond the limit is %g, %.3g orders of magnitude %s %s.",xI(k),dmax,Direction(Sign),LimitName) ;
else
    d=Sign*(xI-LI) ;
    [dmax,k]=max(d) ;
    T=sprintf(" The value furthest beyond the limit is %g, %g %s %s.",xI(k),dmax,Direction(Sign),LimitName) ;
end

end

%%
function T=Direction(Sign)
if Sign<0 ; T="below" ; else ; T="above" ; end
end
