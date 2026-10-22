function Results=TestuvhJacobian(varargin)

%% Taylor test of the Jacobian of the uvh system
%
%   Results=TestuvhJacobian(FileName)                                   % all states saved in FileName (see SaveuvhNewtonState.m)
%   Results=TestuvhJacobian(FileName,Options)
%   Results=TestuvhJacobian(UserVar,CtrlVar,MUA,F0,F1,l1,BCs1)          % a single state
%   Results=TestuvhJacobian(UserVar,CtrlVar,MUA,F0,F1,l1,BCs1,Options)
%
% (15 Oct 2026) For a state x (velocities and thickness at the end of the time step), with residual R(x) and Jacobian K(x) as
% returned by uvhAssembly.m, and a direction d, the remainder
%
%     e(epsilon) = || R(x+epsilon d) - R(x) - epsilon K d ||
%
% is calculated for a sequence of epsilon values (decades). If the Jacobian is consistent with the residual, e decreases as epsilon^2
% (by a factor of 100 per decade), until round-off errors dominate. If a term is missing or wrong, e decreases only as epsilon (by a
% factor of 10 per decade). The test is done separately for the rows of the u, v and h equations, and for several directions:
%
%   "Newton" : the Newton direction of the iteration (if available)
%   "u","v","h" : smooth random perturbations of u, v or h only, scaled to the local magnitude of the field
%   "uvh"    : smooth random perturbations of all three fields
%
% For each direction and block, the observed order is the order between the two smallest epsilon values for which the remainder is
% still clearly above the estimated round-off level. The relative error is e/||epsilon K d|| at that epsilon. If the order is 1, this
% is an estimate of the relative size of the missing or wrong term.
%
% Options (all optional):
%   .Directions    directions to test (default ["Newton","u","v","h","uvh"])
%   .Epsilon       epsilon values, in decades (default 10.^(-1:-1:-9))
%   .States        indices of the saved states to test (file input only; default all)
%   .CtrlVar       structure with fields of CtrlVar to be changed for the test, eg struct("Enforce_bAboveB",false)
%   .NewtonDirection  the Newton direction [du;dv;dh] (single-state input only)
%   .Seed          seed of the random directions (default 1)
%   .MapDirection  direction for which the nodal remainder is kept in Results (default "uvh")
%   .Plot          plot the remainder against epsilon, and the nodal remainder (default false)
%   .Print         print a summary (default true)
%
% Results is a structure array, one element per state, with the context of the state and, for each direction, the remainders, the
% observed orders and the relative errors.
%
% The smallest relative error over all epsilon is the decisive measure: a missing or wrong term of relative size delta gives a relative
% error of at least delta for every epsilon, whereas for a consistent Jacobian it decreases with epsilon until round-off dominates. A low
% order is only flagged if the smallest relative error is also above 1e-4; otherwise the low order is due to round-off noise (eg in the
% h equation for thick ice, where the terms are large).
%
% "lin" in the summary means that the remainder is at round-off level for all epsilon, ie the block is linear in this direction (eg the
% h equation for a perturbation of u or h only, as the flux is bilinear in u and h). This is consistent with the Jacobian.
%
% Note: a perturbation that crosses a non-smooth point of the residual (eg CtrlVar.Enforce_bAboveB) gives first order locally. The
% nodal remainder (Results.Map) shows where the remainder is large.
%
%%

if ischar(varargin{1}) || isstring(varargin{1})
    FileName=string(varargin{1}) ;
    Options=struct() ; if nargin>=2 ; Options=varargin{2} ; end
    Results=TestFile(FileName,Options) ;
else
    narginchk(7,8)
    [UserVar,CtrlVar,MUA,F0,F1,l1,BCs1]=varargin{1:7} ;
    Options=struct() ; if nargin>=8 ; Options=varargin{8} ; end
    dNewton=[] ; if isfield(Options,"NewtonDirection") ; dNewton=Options.NewtonDirection ; end
    Info=struct("State",1,"RunStep",CtrlVar.CurrentRunStepNumber,"time",CtrlVar.time,"Solve",NaN,"Iteration",NaN,"r0",NaN) ;
    if isfield(Options,"Info") ; Info=Options.Info ; end
    Results=TestState(UserVar,CtrlVar,MUA,F0,F1,l1,BCs1,dNewton,Info,Options) ;
    PrintSummary(Results,Options) ;
end

end

%%
function Results=TestFile(FileName,Options)

w=whos("-file",FileName) ; names=string({w.name}) ;
StateNames=sort(names(startsWith(names,"State_"))) ;
if isfield(Options,"States") && ~isempty(Options.States)
    StateNames=StateNames(Options.States) ;
end

MUAname="" ; MUA=[] ; Results=[] ;
for k=1:numel(StateNames)
    L=load(FileName,StateNames(k)) ; S=L.(StateNames(k)) ;
    mn=string(sprintf("MUA_%03d",S.MeshIndex)) ;
    if mn~=MUAname
        M=load(FileName,mn) ; MUA=M.(mn) ; MUAname=mn ;
    end
    Info=struct("State",str2double(extractAfter(StateNames(k),"State_")),"RunStep",S.RunStep,"time",S.time,...
        "Solve",S.Solve,"Iteration",S.Iteration,"r0",S.r0) ;
    R=TestState(S.UserVar,S.CtrlVar,MUA,S.F0,S.F1,S.l1,S.BCs1,S.duvh,Info,Options) ;
    Results=[Results R] ; %#ok<AGROW>
end

PrintSummary(Results,Options) ;

end

%%
function Res=TestState(UserVar,CtrlVar,MUA,F0,F1,l1,BCs1,dNewton,Info,Options)

Directions=["Newton","u","v","h","uvh"] ; if isfield(Options,"Directions") ; Directions=string(Options.Directions) ; end
Epsilon=10.^(-1:-1:-9) ;                if isfield(Options,"Epsilon")    ; Epsilon=Options.Epsilon(:)' ; end
Seed=1 ;                                if isfield(Options,"Seed")       ; Seed=Options.Seed ; end
MapDirection="uvh" ;                    if isfield(Options,"MapDirection") ; MapDirection=string(Options.MapDirection) ; end
DoPlot=false ;                          if isfield(Options,"Plot")       ; DoPlot=Options.Plot ; end

if isfield(Options,"CtrlVar")
    fn=fieldnames(Options.CtrlVar) ;
    for i=1:numel(fn) ; CtrlVar.(fn{i})=Options.CtrlVar.(fn{i}) ; end
end
if isfield(MUA,"dM") ; MUA.dM=[] ; end

RunInfo=UaRunInfo ;
CtrlVar.uvhMatrixAssembly.ZeroFields=false ; CtrlVar.uvhMatrixAssembly.Ronly=false ;
[~,~,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
CtrlVarR=CtrlVar ; CtrlVarR.uvhMatrixAssembly.Ronly=true ;

N=MUA.Nnodes ; rows={1:N,N+1:2*N,2*N+1:3*N} ; BlockNames=["u","v","h"] ;
x=[F1.ub;F1.vb;F1.h] ;

% a priori estimate of the round-off level of the residual in each block (for information; the observed level is used below)
Floor=zeros(3,1) ;
for b=1:3
    Floor(b)=100*eps*(norm(abs(K(rows{b},:))*abs(x))+norm(R(rows{b}))) ;
end

% smooth random fields: random nodal values, averaged over neighbouring nodes five times, scaled to a maximum of one
rng(Seed) ;
con=MUA.connectivity ; nod=size(con,2) ;
I=repmat(con,1,nod) ; J=repelem(con,1,nod) ;
A=spones(sparse(I(:),J(:),1,N,N)) ;
Smooth=spdiags(1./full(sum(A,2)),0,N,N)*A ;
xi=zeros(N,3) ;
for k=1:3
    v=randn(N,1) ; for s=1:5 ; v=Smooth*v ; end ; xi(:,k)=v/max(abs(v)) ;
end
uRef=0.01*max(hypot(F1.ub,F1.vb))+eps ;
du=(abs(F1.ub)+uRef).*xi(:,1) ;
dv=(abs(F1.vb)+uRef).*xi(:,2) ;
dh=max(F1.h,CtrlVar.ThickMin).*xi(:,3) ;
Z=zeros(N,1) ;

Res=Info ;
Res.Floor=Floor ;
Res.Directions=struct("Name",{},"Epsilon",{},"E",{},"Lin",{},"p",{},"Order",{},"RelError",{},"EpsUsed",{},"FloorObs",{},"MinRelError",{}) ;
Res.Map=[] ;

nE=numel(Epsilon) ;
for kd=1:numel(Directions)
    switch Directions(kd)
        case "Newton"
            if isempty(dNewton) ; continue ; end
            % scaled like the random directions (largest relative change one), as the Newton step can be very small
            d=dNewton(:) ;
            sc=max([max(abs(d(rows{1}))./(abs(F1.ub)+uRef)) max(abs(d(rows{2}))./(abs(F1.vb)+uRef)) max(abs(d(rows{3}))./max(F1.h,CtrlVar.ThickMin))]) ;
            if sc>0 ; d=d/sc ; end
        case "u"      ; d=[du;Z;Z] ;
        case "v"      ; d=[Z;dv;Z] ;
        case "h"      ; d=[Z;Z;dh] ;
        case "uvh"    ; d=[du;dv;dh] ;
        otherwise     ; error("TestuvhJacobian:Direction","unknown direction %s",Directions(kd)) ;
    end

    Kd=K*d ;
    E=zeros(3,nE) ; Lin=zeros(3,nE) ; KeepRem=(Directions(kd)==MapDirection) ; RemAll=[] ;
    if KeepRem ; RemAll=zeros(3*N,nE) ; end

    for i=1:nE
        e=Epsilon(i) ;
        F1p=F1 ;
        F1p.ub=F1.ub+e*d(rows{1}) ; F1p.vb=F1.vb+e*d(rows{2}) ; F1p.h=F1.h+e*d(rows{3}) ;
        [~,~,Rp]=uvhAssembly(UserVar,RunInfo,CtrlVarR,MUA,F0,F1p,l1,BCs1) ;
        Rem=Rp-R-e*Kd ;
        for b=1:3
            E(b,i)=norm(Rem(rows{b})) ; Lin(b,i)=norm(e*Kd(rows{b})) ;
        end
        if KeepRem ; RemAll(:,i)=Rem ; end
    end

    % Round-off level: the remainder at the smallest epsilon, where the Taylor terms are negligible (this replaces the a priori
    % estimate, which is too pessimistic for the h equation). Observed order between consecutive epsilon values, where the remainder
    % is clearly above the round-off level.
    FloorObs=max(E(:,end),realmin) ;
    p=nan(3,nE-1) ; Order=nan(3,1) ; RelError=nan(3,1) ; EpsUsed=nan(3,1) ;
    MinRelError=min(E./max(Lin,realmin),[],2) ;   % smallest relative error over all epsilon: a missing term of relative size delta gives at least delta
    for b=1:3
        for i=1:nE-1
            if E(b,i+1)>10*FloorObs(b) && E(b,i)>0
                p(b,i)=log(E(b,i)/E(b,i+1))/log(Epsilon(i)/Epsilon(i+1)) ;
            end
        end
        iv=find(~isnan(p(b,:)),1,"last") ;
        if ~isempty(iv)
            Order(b)=p(b,iv) ; EpsUsed(b)=Epsilon(iv+1) ; RelError(b)=E(b,iv+1)/max(Lin(b,iv+1),realmin) ;
        end
    end

    Res.Directions(end+1)=struct("Name",Directions(kd),"Epsilon",Epsilon,"E",E,"Lin",Lin,"p",p,"Order",Order,...
        "RelError",RelError,"EpsUsed",EpsUsed,"FloorObs",FloorObs,"MinRelError",MinRelError) ;

    if KeepRem
        % nodal remainder, relative to the linear term, at the smallest epsilon above round-off for the u block
        iu=find(Epsilon==EpsUsed(1),1) ; if isempty(iu) ; iu=1 ; end
        Rem=RemAll(:,iu) ; KdE=Epsilon(iu)*Kd ;
        Res.Map=struct("Direction",Directions(kd),"Epsilon",Epsilon(iu),...
            "uv",hypot(Rem(rows{1}),Rem(rows{2}))/max(norm([KdE(rows{1});KdE(rows{2})],Inf),realmin),...
            "h",abs(Rem(rows{3}))/max(norm(KdE(rows{3}),Inf),realmin)) ;
    end
end

if DoPlot
    FigName=sprintf("TestuvhJacobian: state %i",Info.State) ;
    fig=FindOrCreateFigure(FigName) ; clf(fig) ;
    nd=numel(Res.Directions) ;
    for kd=1:nd
        subplot(1,nd,kd) ;
        D=Res.Directions(kd) ;
        loglog(D.Epsilon,D.E','o-') ; hold on
        e1=D.Epsilon(1) ; ref=max(D.E(:,1)) ;
        loglog(D.Epsilon,ref*(D.Epsilon/e1).^2,"k--") ; loglog(D.Epsilon,ref*(D.Epsilon/e1),"k:") ;
        for b=1:3 ; yline(D.FloorObs(b),":") ; end
        xlabel("\epsilon") ; ylabel("remainder") ; title(D.Name) ;
        legend([BlockNames,"\epsilon^2","\epsilon"],Location="southeast") ;
    end
    if ~isempty(Res.Map)
        fig=FindOrCreateFigure(FigName+": nodal remainder") ; clf(fig) ;
        subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,log10(Res.Map.uv+realmin)) ; title(sprintf("log10 |remainder| (u,v), \\epsilon=%g",Res.Map.Epsilon)) ;
        subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,log10(Res.Map.h+realmin)) ; title(sprintf("log10 |remainder| (h), \\epsilon=%g",Res.Map.Epsilon)) ;
    end
end

end

%%
function PrintSummary(Results,Options)

if isfield(Options,"Print") && ~Options.Print ; return ; end

fprintf("\n TestuvhJacobian: observed order of the Taylor remainder (2: consistent, 1: missing or wrong term), and the smallest relative error\n") ;
fprintf(" over all epsilon, ||R(x+eps d)-R(x)-eps K d||/||eps K d|| (a missing term of relative size delta gives at least delta)\n") ;
for k=1:numel(Results)
    R=Results(k) ;
    fprintf(" state %i: run step %i, t=%g, solve %g, Newton iteration %g, r0=%g\n",R.State,R.RunStep,R.time,R.Solve,R.Iteration,R.r0) ;
    for kd=1:numel(R.Directions)
        D=R.Directions(kd) ;
        Flag="" ;
        Low=D.Order<1.5 ;
        if any(Low & D.MinRelError>1e-4)
            Flag="  <-- order below 1.5 and relative error above 1e-4" ;
        elseif any(Low)
            Flag="  (low order only where the remainder is at the noise level; smallest relative error below 1e-4)" ;
        end
        if all(isnan(D.Order)) ; Flag="  (remainder at round-off level for all epsilon, ie linear in this direction)" ; end
        % NaN (remainder at round-off level for all epsilon) is shown as "lin": the block is linear in this direction
        os=compose("%5.2f",D.Order) ; os(isnan(D.Order))="  lin" ;
        ms=compose("%8.1e",D.MinRelError) ;
        fprintf("   %-7s order: u %s  v %s  h %s | smallest rel. error: u %s  v %s  h %s%s\n",D.Name,os,ms,Flag) ;
    end
end

end