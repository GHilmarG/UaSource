function Results=TestInverseHessianEndToEnd(CtrlVar,MUA,F,BCs,l,Meas,Options)

%% End-to-end test of the Hessian of the inverse problem
%
%   Results=TestInverseHessianEndToEnd()
%   Results=TestInverseHessianEndToEnd(CtrlVar,MUA,F,BCs,l)
%   Results=TestInverseHessianEndToEnd(CtrlVar,MUA,F,BCs,l,Meas)
%   Results=TestInverseHessianEndToEnd(...,Name=Value)
%
% Tests that the Hessian calculated by CalcDirectAdjointHessian.m is consistent with the adjoint gradient calculated by
% Misfit.m. For a perturbation $d$ of the inverted fields, $p=(\log_{10} A, B, \log_{10} C)$, the test compares
%
% $$ H \, d  \quad \mathrm{with} \quad \frac{g(p+s\,d)-g(p-s\,d)}{2 s} , $$
%
% where $g=dI/dp$ is the adjoint gradient of the data-misfit function $I$ and $H$ is the direct-adjoint Hessian. Each evaluation
% of $g$ includes a full solution of the forward problem (and of the adjoint problem), and the sequence of calls is the one
% used by JGH.m: p2F.m, uv.m and Misfit.m. The test is therefore independent of every individual term making up the Hessian, and
% it tests all of them together, including the cross terms between the fields and the terms arising from dh/dt measurements.
%
%% What is tested
%
% The test is done for a number of cases, where each case is a combination of inverted fields and measurement types, for
% example inverting for $A$ and $C$ using velocity measurements, or for $A$, $B$ and $C$ using velocity and dh/dt measurements.
%
% For each case the following is done:
%
% * The adjoint gradient is compared with finite differences of the cost function, $g^T d \approx [I(p+s\,d)-I(p-s\,d)]/(2s)$. This is
%   needed because the Hessian is tested against the gradient: a Hessian that is consistently derived from a wrong gradient is
%   consistent with it, and would otherwise pass. The two tests together show that $H$ is the second derivative of the cost function.
%
% * The perturbation directions are pure directions, i.e. $d_A$, $d_B$ and $d_C$ each perturb only one of the fields. The
%   finite-difference estimate of $H d_j$ is compared with the corresponding product from the Hessian, block by block, so that
%   every diagonal and every cross block $(i,j)$ is tested on its own.
%
% * Each comparison is done for a sequence of decreasing steps $s$ (each a factor of two smaller, see StepA, StepB, StepC and
%   NumSteps). A block passes if its relative error at the smallest step is below Tolerance, and if the error has decreased
%   from the largest step. For a correct Hessian the error decreases as $s^2$, so the observed order of convergence (printed in
%   brackets) should be close to 2 until the error reaches the noise floor set by the accuracy of the forward solutions.
%
% * Mutation checks (optional, on by default). Each term of the Hessian is, in turn, left out, and the error is recalculated. A
%   term is "detected" if leaving it out would make the test fail, i.e. if the error without the term exceeds Tolerance. This
%   shows that the test is actually sensitive to that term. A term that is not detected is below the resolution of the test for
%   that case. This is not a failure, but it means that the test says nothing about that term for that case. Making the steps
%   smaller (for B, the closure of the geometry is nonlinear, and the error is then dominated by the step size) and the data
%   weights larger helps.
%
% The relative error of block $i$ is $\| (H d_j)_i - \mathrm{FD}_i \| / \max( \| \mathrm{FD}_i \|, 10^{-3} \| \mathrm{FD} \| )$,
% i.e. blocks that are smaller than a thousandth of the full vector are normalised with that value.
%
%% Which part of the Hessian is tested
%
% The data-misfit part of the Hessian. The regularisation is not included: it is not part of the Misfit.m gradient, and in this test
% the matrices MUA.QA, MUA.QB and MUA.QC are set to zero. This means that the only contribution of Jpp.m is the explicit dependence
% of the dh/dt misfit on $B$ (JhdotBB.m).
%
%% Inputs
%
% With no input arguments a small synthetic problem is generated: a square domain with a floating/grounded transition, a bed
% below sea level, spatially varying A, C and (default) density, no boundary conditions, and synthetic velocity and dh/dt
% measurements obtained from the forward solution plus noise. This is self-contained and can be used as a regression test.
%
% To test the Hessian on a particular problem (mesh, geometry, density, sliding-law coefficients, boundary conditions), supply
% CtrlVar, MUA, F, BCs and l, for example read from a restart file. The measurements are then, by default, synthetic ones created
% from the forward solution at the given parameters. Alternatively, supply Meas, with the fields us, vs, usCov, vsCov, and, if
% dh/dt measurements are used, dhdt and dhdtCov. The covariance matrices of the dh/dt measurements must be diagonal.
%
% The fields F.as and F.ab (if empty) are set to zero. Meas.s, Meas.as and Meas.ab are set from F.
%
% Name-value options:
%
%   Level               "quick" (default) or "full". "quick" runs five cases, "full" runs all combinations of the seven field
%                       combinations with velocity, and with velocity and dh/dt measurements (14 cases).
%   Cases               n x 2 string array, each row being [InvertFor Measurements], e.g. ["-logAGlen-B-" "-uv-dhdt-"].
%                       Overrides Level.
%   Tolerance           Pass tolerance for the relative error of each block (default 1e-4).
%   StepA, StepB, StepC Smallest step sizes $s$. A and C are in units of $\log_{10}$, B in units of the amplitude of the
%                       direction (default 5e-4, 1e-4, 5e-4).
%   NumSteps            Number of steps, each a factor of two larger than the previous one (default 3).
%   AmplitudeA, B, C    Amplitudes of the smooth perturbation directions: 1 in log10(A), 20 m in B, 1 in log10(C).
%   Mutations           Do mutation checks (default true).
%   ErrorsVary          Spatially varying measurement errors (default true). If false, constant errors are used.
%   RhoVaries           Synthetic problem only: spatially varying ice density (default true).
%   dhdtErrorFraction   Synthetic measurements only: dh/dt error as a fraction of the standard deviation of dh/dt (default 0.02).
%   MeshResolution      Synthetic problem only: number of squares in each direction (default 10, giving 221 nodes).
%   Seed                Seed of the random number generator (default 1). The state of the generator is restored on exit.
%   MaxNodes            The Hessian is a full matrix with (number of fields x nodes)^2 entries. The test refuses to run if the
%                       number of nodes is larger than MaxNodes (default 1500).
%   Verbose             Print results (default true).
%
%% Output
%
% Results is a structure with the fields
%
%   Passed      true if all cases passed (for each case both the gradient check and the Hessian check must pass)
%   Cases       structure array with one element for each case, including the relative errors of all blocks and steps
%               (Errors, with dimensions field of output x field of direction x step), the observed orders of convergence and
%               the results of the mutation checks
%   Options, Synthetic, Time
%
%% Nodes with thin ice
%
% The geometrical closure used in B inversions, Calc_bh_From_sBS.m, applies a smooth floor to the thickness, so that h>CtrlVar.ThickMin, and the derivatives with respect
% to B, calculated in dGeometrydB.m, take the floor into account. uv.m then never resets the thickness, and the forward model is a smooth function of B also where the
% ice is thin. No special treatment of thin ice is needed in the test, and all nodes are tested.
%
% If the floor is switched off, CtrlVar.ThickMinWidthRelative=0, this is not so. uv.m then resets thicknesses smaller than ThickMin (and rebuilds b and s at all nodes),
% the forward model is not differentiable with respect to B at nodes with h<=ThickMin, and the derivatives with respect to B (the adjoint gradient and the Hessian) are
% only approximate at such nodes. The test then fails for cases involving B, with errors that do not decrease with the step size, and a note is printed with the number
% of nodes concerned. This is not an error in the Hessian.
%
%% Limitations
%
% The forward-model terms of the Hessian involving the sliding law are implemented for the Weertman sliding law only, and for B,
% CtrlVar.Hh0 must be zero and melange model physics is not supported. The test stops with an error if these conditions are not
% met. Box transformation and the Cholesky mapping of the cost function are not tested.
%
% The test has so far only been run on meshes of three-node elements.
%
% The cost grows quickly with the number of nodes: the sensitivities of the velocities with respect to all inverted fields are
% calculated explicitly. For a case with 1175 nodes, and with A, B and C, one case takes of the order of 20 s, and for the
% "quick" level the test takes about a minute.
%
% Examples:
%
%   TestInverseHessianEndToEnd                                    % self-contained synthetic test
%
%   TestInverseHessianEndToEnd(Level="full")                      % all combinations of fields and measurements
%
%   TestInverseHessianEndToEnd(Cases=["-logAGlen-B-logC-" "-uv-dhdt-"],Tolerance=1e-5,StepB=5e-5)
%
%   TestInverseHessianEndToEnd(CtrlVar,MUA,F,BCs,l)               % on a given problem
%
%  see also: CalcDirectAdjointHessian, Misfit, JGH, Jqq, Jpp, JhdotBB, JhdotBq, Hess_qp, Fqq, Fpp, FAA, FBB, FCC, FAB, FBC
%
%%

arguments
    CtrlVar struct = struct()
    MUA struct = struct()
    F = []
    BCs = []
    l = []
    Meas struct = struct()
    Options.Level (1,1) string {mustBeMember(Options.Level,["quick","full"])} = "quick"
    Options.Cases (:,2) string = strings(0,2)
    Options.Tolerance (1,1) double {mustBePositive} = 1e-4
    Options.StepA (1,1) double {mustBePositive} = 5e-4
    Options.StepB (1,1) double {mustBePositive} = 1e-4
    Options.StepC (1,1) double {mustBePositive} = 5e-4
    Options.NumSteps (1,1) double {mustBeInteger,mustBePositive} = 3
    Options.AmplitudeA (1,1) double = 1
    Options.AmplitudeB (1,1) double = 20
    Options.AmplitudeC (1,1) double = 1
    Options.Mutations (1,1) logical = true
    Options.ErrorsVary (1,1) logical = true
    Options.RhoVaries (1,1) logical = true
    Options.dhdtErrorFraction (1,1) double {mustBePositive} = 0.02
    Options.MeshResolution (1,1) double {mustBeInteger,mustBePositive} = 10
    Options.Seed (1,1) double = 1
    Options.MaxNodes (1,1) double {mustBePositive} = 1500
    Options.Verbose (1,1) logical = true
end

tStart=tic ;
rngState=rng ;
cleanupRng=onCleanup(@() rng(rngState)) ;
rng(Options.Seed) ;

%% The problem

if isempty(fieldnames(CtrlVar))

    Synthetic=true ;
    [CtrlVar,MUA,F,BCs,l]=SyntheticProblem(Options) ;

else

    Synthetic=false ;
    if isempty(fieldnames(MUA)) || isempty(F) || isempty(BCs) || isempty(l)
        error("TestInverseHessianEndToEnd:MissingInput",...
            "CtrlVar, MUA, F, BCs and l must all be supplied. Alternatively, call without any input arguments to use a synthetic test problem.")
    end

end

if SlidingLawCode(CtrlVar.SlidingLaw)~=1   % Weertman
    error("TestInverseHessianEndToEnd:SlidingLaw",...
        "The Hessian is only implemented for the Weertman sliding law, but CtrlVar.SlidingLaw=%s.",string(CtrlVar.SlidingLaw))
end

if CtrlVar.Hh0~=0
    error("TestInverseHessianEndToEnd:Hh0","CtrlVar.Hh0 must be zero, but CtrlVar.Hh0=%g.",CtrlVar.Hh0)
end

N=MUA.Nnodes ;
if N>Options.MaxNodes
    error("TestInverseHessianEndToEnd:TooManyNodes",...
        "The mesh has %i nodes. The Hessian is a full matrix, and the test is limited to %i nodes. Use the option MaxNodes to change this limit.",N,Options.MaxNodes)
end

CtrlVar.InfoLevel=0 ;
CtrlVar.InfoLevelNonLinIt=0 ;

if ~isfield(MUA,"M") || isempty(MUA.M)
    MUA.M=MassMatrix2D1dof(MUA) ;
end

if isempty(F.as) ; F.as=zeros(N,1) ; end
if isempty(F.ab) ; F.ab=zeros(N,1) ; end

% widths of the smooth clipping of A and C: set once, at the start of a run, in Ua2D.m
if ~isfield(CtrlVar,"CminWidth") || ~isfield(CtrlVar,"AGlenminWidth")
    evalc('CtrlVar=SetSmoothFloorWidths(CtrlVar,F) ;') ;
end

% forward solution at the given parameters, which is the starting point of all the other solutions, and which is used to create
% synthetic measurements
evalc('[~,~,F,l]=uv([],[],CtrlVar,MUA,BCs,F,l) ;') ;

if isempty(fieldnames(Meas))
    Meas=SyntheticMeasurements(CtrlVar,MUA,BCs,F,Options) ;
end

Meas.s=F.s ;
Meas.as=F.as ;
Meas.ab=F.ab ;

BCsAdjoint=AdjointBCs(BCs) ;


if isempty(Meas.dhdt)

    N=MUA.Nnodes;
    Meas.dhdt=zeros(N,1) ;
    err=ones(N,1);
    Meas.dhdtCov=spdiags(err.^2,0,N,N) ;
end


%% The cases

Cases=Options.Cases ;

if isempty(Cases)

    switch Options.Level

        case "quick"

            Cases=["-logAGlen-" "-uv-" ;
                "-logC-" "-uv-" ;
                "-logAGlen-logC-" "-uv-" ;
                "-B-" "-uv-dhdt-" ;
                "-logAGlen-B-logC-" "-uv-dhdt-"] ;

        case "full"

            Fields=["-logAGlen-" "-logC-" "-logAGlen-logC-" "-B-" "-logAGlen-B-" "-B-logC-" "-logAGlen-B-logC-"] ;
            Measurements=["-uv-" "-uv-dhdt-"] ;
            Cases=strings(0,2) ;
            for iF=1:numel(Fields)
                for iM=1:numel(Measurements)
                    Cases(end+1,:)=[Fields(iF) Measurements(iM)] ;  %#ok<AGROW>
                end
            end

    end

end

nCases=size(Cases,1) ;

if Options.Verbose
    fprintf("\n TestInverseHessianEndToEnd: %s problem with %i nodes and %i elements, %i case(s), tolerance %g, %i steps.\n",...
        ternary(Synthetic,"synthetic","supplied"),N,MUA.Nele,nCases,Options.Tolerance,Options.NumSteps)
end

CaseResults=cell(1,nCases) ;

for k=1:nCases

    CaseResults{k}=RunOneCase(CtrlVar,MUA,F,BCs,BCsAdjoint,l,Meas,Cases(k,1),Cases(k,2),Options) ;

    if Options.Verbose
        PrintCase(k,nCases,CaseResults{k},Options)
    end

end

Results.Cases=[CaseResults{:}] ;
Results.Passed=all([Results.Cases.Passed]) ;
Results.Options=Options ;
Results.Synthetic=Synthetic ;
Results.Time=toc(tStart) ;

if Options.Verbose

    fprintf("\n %i of %i cases passed (%.0f s). ",nnz([Results.Cases.Passed]),nCases,Results.Time)
    if Results.Passed
        fprintf("The Hessian is consistent with the adjoint gradient in all cases.\n")
    else
        fprintf("FAILED in the following case(s):\n")
        for k=find(~[Results.Cases.Passed])
            fprintf("    invert for %s, measurements %s\n",Results.Cases(k).InvertFor,Results.Cases(k).Measurements)
        end
    end

    if Options.Mutations
        fprintf(" A term whose omission is not detected is below the resolution of the test for that case (see above).\n")
    end

end

if nargout==0
    clear Results
end

end


%% ------------------------------------------------------------------------------------------------------------------------

function Case=RunOneCase(CtrlVar,MUA,F,BCs,BCsAdjoint,l,Meas,InvertFor,Measurements,Options)

tCase=tic ;

CtrlVar.Inverse.InvertFor=InvertFor ;
CtrlVar.Inverse.Measurements=Measurements ;
CtrlVar.Inverse.CalcGradI=true ;
CtrlVar.Inverse.CalcGrad=true ;
CtrlVar.Inverse.DataMisfit.GradientCalculation="adjoint" ;
CtrlVar.Inverse.TestDirectAdjoint.isTrue=false ;
CtrlVar.Inverse.BoxTransform=false ;
CtrlVar.Inverse.CholeskyMappingOfCostFunctionAndGradient=false ;
CtrlVar.Inverse.HessianTerms="-all-" ;
CtrlVar.MapOldToNew.Test=false ;

[isA,isB,isC]=isABC(CtrlVar) ;
[~,isdhdt]=is_uv_dhdt_Meas(CtrlVar) ;

Active=[isA isB isC] ;
N=MUA.Nnodes ;

% nodes with thin ice, see the section on nodes with thin ice in the header: only possible if the floor on the thickness is switched off
NThinIce=0 ;
if isB
    [~,hClosure]=Calc_bh_From_sBS(CtrlVar,MUA,Meas.s,F.B,F.S,F.rho,F.rhow) ;
    NThinIce=nnz(hClosure<=CtrlVar.ThickMin*(1+1e-6)) ;
end

pos=zeros(1,3) ; n=0 ;
for k=1:3
    if Active(k)
        pos(k)=n+1 ;
        n=n+N ;
    end
end
nP=n ;
blk=@(k) pos(k)+(0:N-1) ;     % indices of block k (1=A, 2=B, 3=C) in the vector of inverted fields

% no regularisation: the Hessian of the misfit only
MUA.QA=[] ; MUA.QB=[] ; MUA.QC=[] ; MUA.HBobs=[] ;
Z=sparse(N,N) ;
if isA ; MUA.QA=Z ; end
if isB ; MUA.QB=Z ; end
if isC ; MUA.QC=Z ; end

% parameters and smooth perturbation directions
x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;
xn=(x-min(x))/(max(x)-min(x)) ;
yn=(y-min(y))/(max(y)-min(y)) ;

Pattern{1}=Options.AmplitudeA*sin(2*pi*1.3*xn+0.4).*cos(2*pi*0.9*yn) ;
Pattern{2}=Options.AmplitudeB*sin(2*pi*1.0*xn+1.7).*cos(2*pi*0.7*yn+0.2) ;
Pattern{3}=Options.AmplitudeC*sin(2*pi*0.8*xn).*cos(2*pi*1.1*yn+0.5) ;
Step0=[Options.StepA Options.StepB Options.StepC] ;

p0=zeros(nP,1) ;
if isA ; p0(blk(1))=log10(F.AGlen) ; end
if isB ; p0(blk(2))=F.B ; end
if isC ; p0(blk(3))=log10(F.C) ; end

%% Hessian and adjoint gradient at the base point

[~,g0,Px,Py,Fb,lb]=EvaluateMisfit(CtrlVar,MUA,BCs,BCsAdjoint,l,F,Meas,p0) ;

evalc('H=CalcDirectAdjointHessian(CtrlVar,MUA,Fb,BCs,lb,[],Meas,BCsAdjoint,Px,Py) ;') ;
H=full(H) ;

Symmetry=norm(H-H.',"fro")/norm(H,"fro") ;
IsFinite=all(isfinite(H(:))) ;

%% Comparison with finite differences of the adjoint gradient, for pure directions

nS=Options.NumSteps ;
Errors=nan(3,3,nS) ;        % output block x direction x step (coarsest step first)
Steps=nan(3,nS) ;
GradErrors=nan(3,nS) ;    % adjoint gradient against finite differences of the cost function
Dirs=cell(1,3) ;
FDfine=cell(1,3) ;
Hd=cell(1,3) ;

for j=1:3

    if ~Active(j) ; continue ; end

    d=zeros(nP,1) ;
    d(blk(j))=Pattern{j} ;
    Dirs{j}=d ;
    Hd{j}=H*d ;

    for k=1:nS

        s=Step0(j)*2^(nS-k) ;
        Steps(j,k)=s ;

        [Ip,gp]=EvaluateMisfit(CtrlVar,MUA,BCs,BCsAdjoint,lb,Fb,Meas,p0+s*d) ;
        [Im,gm]=EvaluateMisfit(CtrlVar,MUA,BCs,BCsAdjoint,lb,Fb,Meas,p0-s*d) ;
        FD=(gp-gm)/(2*s) ;

        gd=g0(blk(j)).'*Pattern{j} ;      % directional derivative of the cost according to the adjoint gradient
        GradErrors(j,k)=abs(gd-(Ip-Im)/(2*s))/max(abs(gd),1e-6*norm(g0(blk(j)))*norm(Pattern{j})) ;

        for i=1:3
            if ~Active(i) ; continue ; end
            ref=max(norm(FD(blk(i))),1e-3*norm(FD)) ;
            Errors(i,j,k)=norm(Hd{j}(blk(i))-FD(blk(i)))/ref ;
        end

    end

    FDfine{j}=FD ;      % at the smallest step

end

Tol=Options.Tolerance ;
ErrorFine=Errors(:,:,nS) ;
ErrorCoarse=Errors(:,:,1) ;
Order=log2(ErrorCoarse./ErrorFine)/max(nS-1,1) ;

Pair=(Active.'*Active)>0 ;
BlockPassed=ErrorFine<Tol & (ErrorFine<ErrorCoarse | ErrorCoarse<Tol/10) ;
HessianPassed=all(BlockPassed(Pair)) && IsFinite && Symmetry<1e-8 ;

GradErrorFine=GradErrors(:,nS) ;
GradErrorCoarse=GradErrors(:,1) ;
GradOrder=log2(GradErrorCoarse./GradErrorFine)/max(nS-1,1) ;
GradPassed=all(GradErrorFine(Active)<Tol & (GradErrorFine(Active)<GradErrorCoarse(Active) | GradErrorCoarse(Active)<Tol/10)) ;

Passed=HessianPassed && GradPassed ;

%% Mutation checks

if Options.Mutations
    Mutations=MutationChecks(CtrlVar,MUA,Fb,BCs,BCsAdjoint,lb,Meas,Px,Py,H,Dirs,FDfine,Active,blk,Tol,isdhdt) ;
else
    Mutations=struct("Name",{},"Share",{},"MaxError",{},"Detected",{}) ;
end

Case=struct("InvertFor",InvertFor,"Measurements",Measurements,"Passed",Passed,...
    "Errors",Errors,"ErrorFine",ErrorFine,"ErrorCoarse",ErrorCoarse,"Order",Order,"Steps",Steps,"Active",Active,...
    "GradErrors",GradErrors,"GradErrorFine",GradErrorFine,"GradOrder",GradOrder,"GradPassed",GradPassed,"HessianPassed",HessianPassed,...
    "NThinIce",NThinIce,...
    "Symmetry",Symmetry,"IsFinite",IsFinite,"Mutations",Mutations,"Time",toc(tCase)) ;

end


%% ------------------------------------------------------------------------------------------------------------------------

function Mut=MutationChecks(CtrlVar,MUA,F,BCs,BCsAdjoint,l,Meas,Px,Py,H,Dirs,FDfine,Active,blk,Tol,isdhdt)

%% Leaves out, one at a time, the terms of the Hessian and calculates the resulting error
%
% The terms are first the five terms that CalcDirectAdjointHessian adds up (the two terms involving the sensitivities of the
% velocities, the field-field term, the term from the explicit field dependence of the misfit, and the mixed term). The mixed
% term is obtained as what remains after the other four have been subtracted. For problems involving B further, finer, terms
% are left out: the A-B and B-C cross blocks of the field-field term, and, for dh/dt measurements, the explicit B-velocity term
% of the misfit.

nP=size(H,1) ;
isA=Active(1) ; isB=Active(2) ; isC=Active(3) ;

Hd=cell(1,3) ;
for j=1:3
    if Active(j) ; Hd{j}=H*Dirs{j} ; end
end

Mut=struct("Name",{},"Share",{},"MaxError",{},"Detected",{}) ;
normH=norm(H,"fro") ;

Names=["Jqq   (misfit, second derivatives with respect to u and v)" ;
    "Fqq   (forward model, second derivatives with respect to u and v)" ;
    "Fpp   (forward model, second derivatives with respect to the fields)" ;
    "Jpp   (explicit field dependence of the dh/dt misfit)"] ;
Terms=["-xi Jqq xi-" ; "-xi Fqq xi-" ; "-Fpp-" ; "-Jpp-"] ;

Hrest=H ;

for k=1:numel(Terms)

    Ht=SingleTerm(CtrlVar,MUA,F,BCs,BCsAdjoint,l,Meas,Px,Py,Terms(k)) ;
    Hrest=Hrest-Ht ;

    Share=norm(Ht,"fro")/normH ;
    if Share>0
        MaxError=ErrorWithout(@(d) Ht*d,Hd,Dirs,FDfine,Active,blk) ;
        Mut(end+1)=struct("Name",Names(k),"Share",Share,"MaxError",MaxError,"Detected",MaxError>Tol) ;  %#ok<AGROW>
    end

end

Share=norm(Hrest,"fro")/normH ;
MaxError=ErrorWithout(@(d) Hrest*d,Hd,Dirs,FDfine,Active,blk) ;
Mut(end+1)=struct("Name","Fpq   (mixed field-velocity terms of the forward model and the misfit)","Share",Share,"MaxError",MaxError,"Detected",MaxError>Tol) ;
clear Ht Hrest

%% finer terms involving B

if isB && (isA || isC || isdhdt)

    evalc('[KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC]=duv_hdABC(CtrlVar,MUA,F,l,BCs) ;') ;
    xi=[KdudA KdudB KdudC ; KdvdA KdvdB KdvdC] ;

    if isA

        K=FAB(CtrlVar,MUA,F,BCs,BCsAdjoint,Px,Py) ;
        Apply=@(d) ApplyCross(K,d,blk(1),blk(2),nP) ;
        Mut(end+1)=MutationEntry("Fpp: A-B cross block (FAB)",Apply,Hd,Dirs,FDfine,Active,blk,Tol) ;

    end

    if isC

        K=FBC(CtrlVar,MUA,F,BCs,BCsAdjoint,Px,Py) ;
        Apply=@(d) ApplyCross(K,d,blk(2),blk(3),nP) ;
        Mut(end+1)=MutationEntry("Fpp: B-C cross block (FBC)",Apply,Hd,Dirs,FDfine,Active,blk,Tol) ;

    end

    if isdhdt

        KJq=JhdotBq(CtrlVar,MUA,F,Meas)*xi ;
        Apply=@(d) ApplyRowBlock(KJq,d,blk(2),nP) ;
        Mut(end+1)=MutationEntry("Fpq: explicit B-velocity term of the dh/dt misfit (JhdotBq)",Apply,Hd,Dirs,FDfine,Active,blk,Tol) ;

    end

end

end


function Entry=MutationEntry(Name,Apply,Hd,Dirs,FDfine,Active,blk,Tol)

MaxError=ErrorWithout(Apply,Hd,Dirs,FDfine,Active,blk) ;
% Share is NaN: the size of these finer terms is not comparable with the norms of the full terms
Entry=struct("Name",Name,"Share",NaN,"MaxError",MaxError,"Detected",MaxError>Tol) ;

end


function MaxError=ErrorWithout(Apply,Hd,Dirs,FDfine,Active,blk)

% maximum, over all blocks, of the relative error of H*d if the term Apply(d) is left out
MaxError=0 ;

for j=1:3

    if ~Active(j) ; continue ; end

    r=Hd{j}-Apply(Dirs{j}) ;
    nFD=norm(FDfine{j}) ;

    for i=1:3
        if ~Active(i) ; continue ; end
        ref=max(norm(FDfine{j}(blk(i))),1e-3*nFD) ;
        MaxError=max(MaxError,norm(r(blk(i))-FDfine{j}(blk(i)))/ref) ;
    end

end

end


function Ed=ApplyCross(K,d,rowsA,rowsB,nP)

% product of d with a symmetric matrix having K in block (A,B) and K.' in block (B,A)
Ed=zeros(nP,1) ;
Ed(rowsA)=K*d(rowsB) ;
Ed(rowsB)=K.'*d(rowsA) ;

end


function Ed=ApplyRowBlock(K,d,rows,nP)

% product of d with a symmetric matrix of the form E+E.', where E has the (N x nP) matrix K in the rows 'rows'
Ed=zeros(nP,1) ;
Ed(rows)=K*d ;
Ed=Ed+K.'*d(rows) ;

end


function Ht=SingleTerm(CtrlVar,MUA,F,BCs,BCsAdjoint,l,Meas,Px,Py,Term) %#ok<INUSD>

CtrlVar.Inverse.HessianTerms=Term ;
evalc('Ht=CalcDirectAdjointHessian(CtrlVar,MUA,F,BCs,l,[],Meas,BCsAdjoint,Px,Py) ;') ;
Ht=full(Ht) ;

end


%% ------------------------------------------------------------------------------------------------------------------------

function [I,g,Psi_x,Psi_y,F,l]=EvaluateMisfit(CtrlVar,MUA,BCs,BCsAdjoint,l,F,Meas,p) %#ok<INUSL>

%% Misfit and its adjoint gradient at the parameter vector p
%
% Same sequence as in JGH.m: map the parameters to fields, solve the forward problem, and then calculate the misfit and the
% adjoint gradient. The forward problem is started from the F supplied, which should be a converged solution.

evalc('F=p2F(CtrlVar,MUA,p,F,Meas,[]) ;') ;
evalc('[~,~,F,l,dFduv]=uv([],[],CtrlVar,MUA,BCs,F,l) ;') ;
evalc('[I,g,Psi_x,Psi_y,F,~]=Misfit(CtrlVar,MUA,BCs,F,l,[],Meas,BCsAdjoint,dFduv) ;') ;
g=full(g) ;

end


function BCsAdjoint=AdjointBCs(BCs)

% boundary conditions for the adjoint problem are homogeneous
BCsAdjoint=BCs ;
Props=string(properties(BCs)) ;
for k=1:numel(Props)
    if endsWith(Props(k),"Value") && isnumeric(BCsAdjoint.(Props(k))) && ~isempty(BCsAdjoint.(Props(k)))
        BCsAdjoint.(Props(k))=0*BCsAdjoint.(Props(k)) ;
    end
end

end


function Meas=SyntheticMeasurements(CtrlVar,MUA,BCs,F,Options)

%% Synthetic measurements: the forward solution plus noise
%
% The noise is a fraction (0.3) of the measurement error, which is, if ErrorsVary is true, proportional to the speed and varies
% in space.

N=MUA.Nnodes ;
x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;
xn=(x-min(x))/(max(x)-min(x)) ;
yn=(y-min(y))/(max(y)-min(y)) ;

Speed=hypot(F.ub,F.vb) ;

if Options.ErrorsVary
    SpeedRef=max(Speed,0.5*median(Speed)) ;
    su=0.2*SpeedRef.*(1+0.5*sin(2*pi*xn)) ;
    sv=0.2*SpeedRef.*(1+0.5*cos(2*pi*yn)) ;
else
    su=zeros(N,1)+0.2*median(Speed) ;
    sv=su ;
end

Meas.us=F.ub+0.3*su.*randn(N,1) ;
Meas.vs=F.vb+0.3*sv.*randn(N,1) ;
Meas.usCov=spdiags(su.^2,0,N,N) ;
Meas.vsCov=spdiags(sv.^2,0,N,N) ;

evalc('[~,dhdt]=dhdtExplicit([],CtrlVar,MUA,F,BCs) ;') ;
sh=Options.dhdtErrorFraction*max(std(dhdt),1e-3)*(1+0.5*Options.ErrorsVary*sin(2*pi*xn)) ;
Meas.dhdt=dhdt+0.3*sh.*randn(N,1) ;
Meas.dhdtCov=spdiags(sh.^2,0,N,N) ;

end


function [CtrlVar,MUA,F,BCs,l]=SyntheticProblem(Options)

%% A small synthetic problem
%
% A square domain, 100 km x 100 km, with a bed below sea level and an ice thickness that is within a few tens of metres of
% flotation, so that the domain consists of a mixture of grounded and floating ice, and the grounding line is crossed by many
% elements. A and C vary in space, and, optionally, so does the ice density. There are no boundary conditions: the basal
% friction in the grounded parts makes the problem well-posed.

evalc('CtrlVar=Ua2D_DefaultParameters ;') ;

CtrlVar.TriNodes=3 ;
CtrlVar.MeshGenerator="UaSquareMesh" ;
CtrlVar.UaSquareMesh.xmin=-50e3 ;
CtrlVar.UaSquareMesh.xmax=50e3 ;
CtrlVar.UaSquareMesh.ymin=-50e3 ;
CtrlVar.UaSquareMesh.ymax=50e3 ;
CtrlVar.UaSquareMesh.nx=Options.MeshResolution ;
CtrlVar.UaSquareMesh.ny=Options.MeshResolution ;
CtrlVar.UaSquareMesh.Refine=true ;
CtrlVar.CalcMUA_Derivatives=true ;
CtrlVar.MUA.MassMatrix=true ;
CtrlVar.SlidingLaw="Weertman" ;
CtrlVar.kH=1 ;
CtrlVar.Hh0=0 ;
CtrlVar.GroundingFloatingMaskContains="GF nodes" ;
CtrlVar.InfoLevel=0 ;
CtrlVar.InfoLevelNonLinIt=0 ;

evalc('[coordinates,connectivity]=UaSquareMesh(CtrlVar) ; MUA=CreateMUA(CtrlVar,connectivity,coordinates) ;') ;

N=MUA.Nnodes ;
x=MUA.coordinates(:,1) ;
y=MUA.coordinates(:,2) ;

F=UaFields ;
F.x=x ; F.y=y ;
F.time=0 ; F.dt=0.1 ;
F.g=0.00981 ; F.rhow=1027 ; F.alpha=0 ;
F.rho=920+zeros(N,1) ;
if Options.RhoVaries
    F.rho=920+40*sin(2*pi*x/1.2e5).*cos(2*pi*y/9e4) ;
end
F.S=zeros(N,1) ;
F.n=3+zeros(N,1) ;
F.m=3+zeros(N,1) ;

F.B=-350+250*sin(2*pi*x/8e4+1).*cos(2*pi*y/6e4) ;
hf=F.rhow*(F.S-F.B)./F.rho ;
F.h=hf+30*sin(2*pi*x/4.4e4).*cos(2*pi*y/5.2e4+0.5) ;

[F.b,F.s]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow) ;
[F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow) ;

F.AGlen=2e-9*(1+0.3*sin(2*pi*x/6.8e4)) ;
F.C=3e-3*(1+0.3*cos(2*pi*y/7.6e4+0.3)) ;
F.AGlenmin=1e-20 ; F.AGlenmax=1e20 ;
F.Cmin=1e-20 ; F.Cmax=1e50 ;

F.ub=zeros(N,1)+50 ; F.vb=zeros(N,1) ;
F.ud=zeros(N,1) ; F.vd=zeros(N,1) ;
F.as=zeros(N,1) ; F.ab=zeros(N,1) ;
F.dhdt=zeros(N,1) ;
F.dubdt=zeros(N,1) ; F.dvbdt=zeros(N,1) ; F.duddt=zeros(N,1) ; F.dvddt=zeros(N,1) ;

BCs=BoundaryConditions ;
l=UaLagrangeVariables ;

end


%% ------------------------------------------------------------------------------------------------------------------------

function PrintCase(k,nCases,Case,Options)

Names=["A" "B" "C"] ;

Status="PASSED" ;
if ~Case.Passed
    Status="FAILED:" ;
    if ~Case.GradPassed ; Status=Status+" gradient" ; end
    if ~Case.HessianPassed ; Status=Status+" Hessian" ; end
end

fprintf("\n Case %i/%i:  invert for %s,  measurements %s:  %s  (%.1f s)\n",k,nCases,Case.InvertFor,Case.Measurements,Status,Case.Time)

if Case.NThinIce>0
    fprintf("   Note: %i nodes have an ice thickness at or below ThickMin in the geometrical closure, which is only possible if the floor on the thickness is switched off (CtrlVar.ThickMinWidthRelative=0). The forward model is then not differentiable with respect to B there, and the failure of this case is expected.\n",Case.NThinIce)
end
fprintf("   Relative error of (H d_j)_i against finite differences of the adjoint gradient at the smallest step,\n")
fprintf("   with the observed order of convergence in brackets (expected: 2).\n")

fprintf("   %-10s","")
for j=1:3
    if Case.Active(j) ; fprintf("%-18s","d_"+Names(j)) ; end
end
fprintf("\n")

for i=1:3

    if ~Case.Active(i) ; continue ; end

    fprintf("   output %-3s",Names(i))
    for j=1:3
        if Case.Active(j)
            fprintf("%9.1e (%4.1f)  ",Case.ErrorFine(i,j),Case.Order(i,j))
        end
    end
    fprintf("\n")

end

fprintf("   Adjoint gradient against finite differences of the cost function (g'*d_j against (I(p+s d_j)-I(p-s d_j))/2s):\n   ")
for j=1:3
    if Case.Active(j) ; fprintf("d_%s %9.1e (%4.1f)   ",Names(j),Case.GradErrorFine(j),Case.GradOrder(j)) ; end
end
fprintf("\n")

fprintf("   symmetry of H: %.1e,  finite: %i\n",Case.Symmetry,Case.IsFinite)

if ~isempty(Case.Mutations)

    fprintf("   Mutation checks: the error if a term is left out (the term is detected if the error exceeds the tolerance %g).\n",Options.Tolerance)
    fprintf("   The first column is the norm of the term divided by the norm of H (the terms partly cancel, so these need not add up to one).\n")
    for m=1:numel(Case.Mutations)
        Mut=Case.Mutations(m) ;
        if Mut.Detected
            Det="detected" ;
        else
            Det="below resolution" ;
        end
        if isnan(Mut.Share)
            ShareString="       " ;
        else
            ShareString=sprintf("%6.3f ",Mut.Share) ;
        end
        fprintf("      %s %-72s error without term %9.2e  %s\n",ShareString,Mut.Name,Mut.MaxError,Det)
    end

end

end


function s=ternary(condition,a,b)

if condition
    s=a ;
else
    s=b ;
end

end
