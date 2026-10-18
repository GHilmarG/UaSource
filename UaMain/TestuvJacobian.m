function Results=TestuvJacobian(varargin)

%% Taylor test of the Jacobian of the uv (SSA/SSTREAM) system
%
%   Results=TestuvJacobian(RestartFileName)                 % the state in a restart file (F, MUA, BCs, CtrlVarInRestartFile)
%   Results=TestuvJacobian(RestartFileName,Options)
%   Results=TestuvJacobian(CtrlVar,MUA,F,BCs)               % a given state
%   Results=TestuvJacobian(CtrlVar,MUA,F,BCs,Options)
%
% (17 Oct 2026) The uv counterpart of TestuvhJacobian.m. For a state x=[ub;vb] (with all other fields in F fixed), with residual R(x)
% and Jacobian K(x) as returned by uvMatrixAssembly.m (the same calls as in SSTREAM2dNR2.m), and a direction d, the remainder
%
%     e(epsilon) = || R(x+epsilon d) - R(x) - epsilon K d ||
%
% is calculated for a sequence of epsilon values (decades). If the Jacobian is consistent with the residual, e decreases as epsilon^2
% (by a factor of 100 per decade), until round-off errors dominate. If a term is missing or wrong, e decreases only as epsilon (by a
% factor of 10 per decade). The test is done separately for the rows of the u and v equations, and for several directions:
%
%   "Newton" : the Newton direction at the state, including the boundary conditions (or Options.NewtonDirection)
%   "u","v"  : smooth random perturbations of u or v only, scaled to the local magnitude of the velocity
%   "uv"     : smooth random perturbations of both
%
% For each direction and block, the observed order is the order between the two smallest epsilon values for which the remainder is
% still clearly above the round-off level, and the smallest relative error over all epsilon, e/||epsilon K d||, is the decisive
% measure: a missing or wrong term of relative size delta gives a relative error of at least delta for every epsilon, whereas for a
% consistent Jacobian it decreases with epsilon until round-off dominates. A low order is only flagged if the smallest relative error
% is also above 1e-4.
%
% In addition, the symmetry of the Jacobian is reported, ||K-K'||/||K|| (Frobenius norm): with Glen's flow law and Weertman sliding
% the Newton Jacobian of the SSA is symmetric, and SSTREAM2dNR2.m then uses solveKApeSymmetric.m. The sliding law and the exponents
% n and m are reported as well. (uvMatrixAssemblySSTREAM.m always forms the full Newton Jacobian.)
%
% Consistency is a property of the residual and the Jacobian at any state, so the uv iterates of a solve need not be saved: other
% states can be tested by changing F.ub and F.vb, eg scaling them or adding perturbations, and calling this function again.
%
% Options (all optional):
%   .Directions       directions to test (default ["Newton","u","v","uv"])
%   .Epsilon          epsilon values, in decades (default 10.^(-1:-1:-9))
%   .CtrlVar          structure with fields of CtrlVar to be changed for the test, eg struct("SlidingLaw","Budd")
%   .NewtonDirection  the Newton direction [du;dv] (default: calculated at the state)
%   .Seed             seed of the random directions (default 1)
%   .MapDirection     direction for which the nodal remainder is kept in Results (default "uv")
%   .Label            text describing the state, shown in the summary (default "")
%   .ModifyJacobian   function handle applied to the Jacobian before the test, K=Options.ModifyJacobian(K), to check the test itself,
%                     eg @(K) K+1e-3*spdiags(diag(K),0,size(K,1),size(K,1)) must give first order and a relative error of about 1e-3
%   .Plot             plot the remainder against epsilon, and the nodal remainder (default false)
%   .Print            print a summary (default true)
%
% Results is a structure with the label of the state, the symmetry of the Jacobian, the options that determine the Jacobian, and, for
% each direction, the remainders, the observed orders and the relative errors.
%
% Note: a perturbation that crosses a non-smooth point of the residual (eg a velocity or a strain rate passing through a
% regularisation limit) gives first order locally. The nodal remainder (Results.Map) shows where the remainder is large.
%
%%

if ischar(varargin{1}) || isstring(varargin{1})
    FileName=string(varargin{1}) ;
    Options=struct() ; if nargin>=2 ; Options=varargin{2} ; end
    w0=warning("off","all") ;     % eg the warning that a decomposition object can not be loaded
    S=load(FileName,"F","MUA","BCs","CtrlVarInRestartFile") ;
    warning(w0) ;
    CtrlVar=S.CtrlVarInRestartFile ;
    CtrlVar.Parallel.uvAssembly.spmd.isOn=false ; CtrlVar.Parallel.uvhAssembly.spmd.isOn=false ;
    MUA=UpdateMUA(CtrlVar,S.MUA) ;
    F=S.F ; BCs=S.BCs ;
    if ~isfield(Options,"Label") ; Options.Label=FileName ; end
else
    narginchk(4,5)
    [CtrlVar,MUA,F,BCs]=varargin{1:4} ;
    Options=struct() ; if nargin>=5 ; Options=varargin{5} ; end
end

Results=TestState(CtrlVar,MUA,F,BCs,Options) ;
PrintSummary(Results,Options) ;

end

%%
function Res=TestState(CtrlVar,MUA,F,BCs,Options)

Directions=["Newton","u","v","uv"] ; if isfield(Options,"Directions") ; Directions=string(Options.Directions) ; end
Epsilon=10.^(-1:-1:-9) ;             if isfield(Options,"Epsilon")    ; Epsilon=Options.Epsilon(:)' ; end
Seed=1 ;                             if isfield(Options,"Seed")       ; Seed=Options.Seed ; end
MapDirection="uv" ;                  if isfield(Options,"MapDirection") ; MapDirection=string(Options.MapDirection) ; end
DoPlot=false ;                       if isfield(Options,"Plot")       ; DoPlot=Options.Plot ; end
Label="" ;                           if isfield(Options,"Label")      ; Label=string(Options.Label) ; end

if isfield(Options,"CtrlVar")
    fn=fieldnames(Options.CtrlVar) ;
    for i=1:numel(fn) ; CtrlVar.(fn{i})=Options.CtrlVar.(fn{i}) ; end
end
if isfield(MUA,"dM") ; MUA.dM=[] ; end

% residual and Jacobian, as in SSTREAM2dNR2.m
RunInfo=UaRunInfo ;
CtrlVar.uvAssembly.ZeroFields=false ; CtrlVar.uvMatrixAssembly.Ronly=false ;
[RunInfo,R,K]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F,BCs) ;
CtrlVarR=CtrlVar ; CtrlVarR.uvMatrixAssembly.Ronly=true ;
Modified=isfield(Options,"ModifyJacobian") && ~isempty(Options.ModifyJacobian) ;
if Modified ; K=Options.ModifyJacobian(K) ; end   % only for checking the test itself

N=MUA.Nnodes ; rows={1:N,N+1:2*N} ; BlockNames=["u","v"] ;
x=[F.ub;F.vb] ;

Res=struct() ;
Res.Label=Label ;
Res.Symmetry=norm(K-K',"fro")/max(norm(K,"fro"),realmin) ;
Res.JacobianModified=Modified ;
Res.SlidingLaw="unknown" ; if isfield(CtrlVar,"SlidingLaw") ; Res.SlidingLaw=string(CtrlVar.SlidingLaw) ; end
Res.n=[min(F.n) max(F.n)] ; Res.m=[min(F.m) max(F.m)] ;
Res.ResidualNorm=norm(R) ;

% Newton direction at the state, with the boundary conditions, as in SSTREAM2dNR2.m (the multipliers only shift the multiplier
% increment, so they are not needed for the velocity increment)
dNewton=[] ;
if isfield(Options,"NewtonDirection") && ~isempty(Options.NewtonDirection)
    dNewton=Options.NewtonDirection(:) ;
elseif any(Directions=="Newton")
    MLC=BCs2MLC(CtrlVar,MUA,BCs) ; L=MLC.ubvbL ; cuv=MLC.ubvbRhs ;
    if isempty(L)
        frhs=-R ; grhs=[] ;
    else
        frhs=-R ; grhs=cuv-L*x ;
    end
    if issymmetric(K)
        dNewton=solveKApeSymmetric(K,L,frhs,grhs,zeros(2*N,1),zeros(numel(grhs),1),CtrlVar) ;
    else
        dNewton=solveKApe(K,L,frhs,grhs,zeros(2*N,1),zeros(numel(grhs),1),CtrlVar) ;
    end
end

% smooth random fields: random nodal values, averaged over neighbouring nodes five times, scaled to a maximum of one
rng(Seed) ;
con=MUA.connectivity ; nod=size(con,2) ;
I=repmat(con,1,nod) ; J=repelem(con,1,nod) ;
A=spones(sparse(I(:),J(:),1,N,N)) ;
Smooth=spdiags(1./full(sum(A,2)),0,N,N)*A ;
xi=zeros(N,2) ;
for k=1:2
    v=randn(N,1) ; for s=1:5 ; v=Smooth*v ; end ; xi(:,k)=v/max(abs(v)) ;
end
uRef=0.01*max(hypot(F.ub,F.vb))+eps ;
du=(abs(F.ub)+uRef).*xi(:,1) ;
dv=(abs(F.vb)+uRef).*xi(:,2) ;
Z=zeros(N,1) ;

Res.Directions=struct("Name",{},"Epsilon",{},"E",{},"Lin",{},"p",{},"Order",{},"RelError",{},"EpsUsed",{},"FloorObs",{},"MinRelError",{}) ;
Res.Map=[] ;

nE=numel(Epsilon) ;
for kd=1:numel(Directions)
    switch Directions(kd)
        case "Newton"
            if isempty(dNewton) ; continue ; end
            % scaled like the random directions (largest relative change one), as the Newton step can be very small
            d=dNewton ;
            sc=max([max(abs(d(rows{1}))./(abs(F.ub)+uRef)) max(abs(d(rows{2}))./(abs(F.vb)+uRef))]) ;
            if sc>0 ; d=d/sc ; end
        case "u"      ; d=[du;Z] ;
        case "v"      ; d=[Z;dv] ;
        case "uv"     ; d=[du;dv] ;
        otherwise     ; error("TestuvJacobian:Direction","unknown direction %s",Directions(kd)) ;
    end

    Kd=K*d ;
    E=zeros(2,nE) ; Lin=zeros(2,nE) ; KeepRem=(Directions(kd)==MapDirection) ; RemAll=[] ;
    if KeepRem ; RemAll=zeros(2*N,nE) ; end

    for i=1:nE
        e=Epsilon(i) ;
        Fp=F ;
        Fp.ub=F.ub+e*d(rows{1}) ; Fp.vb=F.vb+e*d(rows{2}) ;
        [~,Rp]=uvMatrixAssembly(RunInfo,CtrlVarR,MUA,Fp,BCs) ;
        Rem=Rp-R-e*Kd ;
        for b=1:2
            E(b,i)=norm(Rem(rows{b})) ; Lin(b,i)=norm(e*Kd(rows{b})) ;
        end
        if KeepRem ; RemAll(:,i)=Rem ; end
    end

    % round-off level: the remainder at the smallest epsilon; observed order between consecutive epsilon values, where the remainder
    % is clearly above the round-off level
    FloorObs=max(E(:,end),realmin) ;
    p=nan(2,nE-1) ; Order=nan(2,1) ; RelError=nan(2,1) ; EpsUsed=nan(2,1) ;
    MinRelError=min(E./max(Lin,realmin),[],2) ;
    for b=1:2
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

    Res.Directions(end+1)=struct("Name",Directions(kd),"Epsilon",Epsilon,"E",E,"Lin",Lin,"p",p,"Order",Order,...   %#ok<AGROW>
        "RelError",RelError,"EpsUsed",EpsUsed,"FloorObs",FloorObs,"MinRelError",MinRelError) ;

    if KeepRem
        % nodal remainder, relative to the linear term, at the smallest epsilon above round-off for the u block
        iu=find(Epsilon==EpsUsed(1),1) ; if isempty(iu) ; iu=1 ; end
        Rem=RemAll(:,iu) ; KdE=Epsilon(iu)*Kd ;
        Res.Map=struct("Direction",Directions(kd),"Epsilon",Epsilon(iu),...
            "uv",hypot(Rem(rows{1}),Rem(rows{2}))/max(norm(hypot(KdE(rows{1}),KdE(rows{2})),Inf),realmin)) ;
    end
end

if DoPlot
    FigName="TestuvJacobian" ; if Label~="" ; FigName=FigName+": "+Label ; end
    fig=FindOrCreateFigure(FigName) ; clf(fig) ;
    nd=numel(Res.Directions) ;
    for kd=1:nd
        subplot(1,nd,kd) ;
        D=Res.Directions(kd) ;
        loglog(D.Epsilon,D.E','o-') ; hold on
        e1=D.Epsilon(1) ; ref=max(D.E(:,1)) ;
        loglog(D.Epsilon,ref*(D.Epsilon/e1).^2,"k--") ; loglog(D.Epsilon,ref*(D.Epsilon/e1),"k:") ;
        for b=1:2 ; yline(D.FloorObs(b),":") ; end
        xlabel("\epsilon") ; ylabel("remainder") ; title(D.Name) ;
        legend([BlockNames,"\epsilon^2","\epsilon"],Location="southeast") ;
    end
    if ~isempty(Res.Map)
        fig=FindOrCreateFigure(FigName+": nodal remainder") ; clf(fig) ;
        PlotMeshScalarVariable(CtrlVar,MUA,log10(Res.Map.uv+realmin)) ;
        title(sprintf("log10 |remainder| (u,v), \\epsilon=%g",Res.Map.Epsilon)) ;
    end
end

end

%%
function PrintSummary(R,Options)

if isfield(Options,"Print") && ~Options.Print ; return ; end

fprintf("\n TestuvJacobian: observed order of the Taylor remainder (2: consistent, 1: missing or wrong term), and the smallest relative error\n") ;
fprintf(" over all epsilon, ||R(x+eps d)-R(x)-eps K d||/||eps K d|| (a missing term of relative size delta gives at least delta)\n") ;
if R.Label~="" ; fprintf(" state: %s\n",R.Label) ; end
fprintf(" sliding law %s, n in [%g %g], m in [%g %g] | ||R||=%g | symmetry of K: ||K-K''||/||K||=%.1e\n",...
    R.SlidingLaw,R.n,R.m,R.ResidualNorm,R.Symmetry) ;
if R.JacobianModified ; fprintf(" NOTE: the Jacobian was modified for the test (Options.ModifyJacobian)\n") ; end
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
    os=compose("%5.2f",D.Order) ; os(isnan(D.Order))="  lin" ;
    ms=compose("%8.1e",D.MinRelError) ;
    fprintf("   %-7s order: u %s  v %s | smallest rel. error: u %s  v %s%s\n",D.Name,os,ms,Flag) ;
end

end
