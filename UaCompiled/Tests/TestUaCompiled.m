classdef TestUaCompiled < matlab.unittest.TestCase

    %% Regression test of the compiled uvh assembly (UaSource/UaCompiled) against the vectorised assembly
    %
    %   results=runtests("TestUaCompiled")                                   % all tests, about 1 minute
    %   table(results)                                                       % one row per test
    %   results=runtests("TestUaCompiled","ProcedureName","testSlidingLaws")  % one group of tests
    %
    % What it tests, on a small synthetic state built here (no files needed: grounded and floating ice with a grounding line, a
    % patch where the thickness goes down to about ThickMin, a level-set field crossing the domain, melange fields, and
    % velocities that are not a solution, so that the residuals are not small remainders of cancelling terms):
    %
    %   - for each supported option (sliding laws with and without melange, SUPG tau options, thickness penalty functions, thickness
    %     barrier, active-set constraints, level-set and mass-balance geometry feedback, ResetThicknessToMinThickness,
    %     Enforce_bAboveB), element type (linear, quadratic, cubic) and quadrature rule (1 to 33 points), that
    %       * the compiled path is actually used (CompiledAssemblyCheck.m returns ok=true; a silent fallback to the vectorised
    %         assembly would otherwise make every comparison pass), and
    %       * the residual (momentum and thickness rows separately) and the Jacobian agree with the vectorised assembly to
    %         rounding errors (Tolerance, relative to the largest value), for the full and the residual-only assembly;
    %   - the basal drag, effective viscosity and floating mask at the integration points (OnlyCalcBasalDragAndEffectiveViscosity);
    %   - that a quadrature rule outside the compiled range (42 points) falls back to the vectorised assembly, with identical results;
    %   - one implicit uvh time step (uvh.m, Newton-Raphson): the same number of iterations, and the same solution to rounding errors;
    %   - that the MEX files are up to date (UaCompiledBuildInfo.m).
    %
    % What it does not test: changes to the physics made in the same way in both assemblies (the test checks that the two agree,
    % not that the physics is right).
    %
    % Run it after any change to uvhElementKernel.m or the other files in UaCompiled, to the functions the kernel compiles (see
    % UaCompiledBuildInfo.m, eg BasalDrag.m, SUPGtau.m), or to the vectorised uvh assembly (uvhMatrixAssembly.m,
    % uvhAssemblyIntPointImplicitSUPG_v2.m, uvhNodalLoopSSTREAM_v2.m), and after rebuilding the MEX files.
    %
    % See also: README.md in UaCompiled, BuildUaCompiled, CompiledAssemblyCheck
    %
    % (26 Oct 2026)
    %%

    properties (Constant)
        Tolerance = 1e-11     % assembly: relative to the largest value of R (momentum and thickness rows separately) and of K
        ToleranceSolve = 1e-9 % one uvh time step: relative difference of h, ub, vb
    end

    properties (TestParameter)
        SlidingLaw = {"Weertman","Budd","Coulomb","Tsai","Cornford","Umbi","Joughin"}
        Melange = struct("withoutMelange",false,"withMelange",true)
        TauOption = {"tau1","tau2","taus","taut","tau0"}
        PenaltyFunction = {"softplus","polynomial","exponential"}
        Option = {"LevelSetFeedback","ThicknessBarrier","ActiveSetConstraints","ResetThicknessToMinThickness","Enforce_bAboveB","MassBalanceGeometryFeedback"}
        TriNodes = struct("linear",3,"quadratic",6,"cubic",10)
        QuadratureDegree = struct("degree1",1,"degree2",2,"degree4",4,"degree8",8,"degree12",12)
    end

    properties
        Base   % the synthetic state on the linear mesh (built once, in TestClassSetup)
    end

    methods (TestClassSetup)

        function setupStateAndCheckMexFiles(testCase)
            for m=["uvhKernelFull_mex","uvhKernelRonly_mex","uvhKernelFull_n3i6_mex","uvhKernelRonly_n3i6_mex","SparseFromCSC"]
                testCase.assertEqual(exist(m,"file"),3, ...
                    "The MEX file "+m+" was not found. Build the MEX files with BuildUaCompiled.m (see README.md in UaCompiled).")
            end
            clear CompiledAssemblyCheck   % start from a fresh state (eg after an earlier failure)
            testCase.Base=SyntheticUvhState(3,40,[]) ;
        end

    end

    methods (Test)

        function testBuildInformation(testCase)
            [ok,Reason,Note]=UaCompiledBuildInfo("check",["uvhKernelFull_mex","uvhKernelRonly_mex","uvhKernelFull_n3i6_mex","uvhKernelRonly_n3i6_mex","SparseFromCSC"]) ;
            testCase.verifyTrue(ok,"MEX files not up to date: "+Reason)
            testCase.verifyEqual(Note,"","No build information: "+Note)
        end

        function testSlidingLaws(testCase,SlidingLaw,Melange)
            S=testCase.Base ;
            S.CtrlVar.SlidingLaw=SlidingLaw ;
            S.CtrlVar.IncludeMelangeModelPhysics=Melange ;
            testCase.verifySameAssembly(S,"sliding law "+SlidingLaw+", melange "+Melange)
        end

        function testTauOptions(testCase,TauOption)
            S=testCase.Base ;
            S.CtrlVar.uvh.SUPG.tau=TauOption ;
            testCase.verifySameAssembly(S,"SUPG tau option "+TauOption)
        end

        function testThicknessPenalty(testCase,PenaltyFunction)
            S=testCase.Base ;
            S.CtrlVar.ThicknessPenalty=true ;
            S.CtrlVar.ThicknessPenaltyNodal=false ;
            S.CtrlVar.ThicknessPenaltyMassBalanceFeedbackFunction=PenaltyFunction ;
            testCase.verifySameAssembly(S,"thickness penalty "+PenaltyFunction)
        end

        function testOptions(testCase,Option)
            S=testCase.Base ;
            switch Option
                case "LevelSetFeedback"
                    S.CtrlVar.LevelSetMethod=true ;
                    S.CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback=true ;
                    S.F1.LSFMask=[] ; S.F0.LSFMask=[] ;   % computed from LSF by uvhMatrixAssembly
                case "ThicknessBarrier"
                    S.CtrlVar.ThicknessBarrier=true ;
                    S.CtrlVar.ThicknessBarrierMassBalanceFeedbackCoeffLog=1 ;
                case "ActiveSetConstraints"
                    I=find(S.F1.h<2*S.CtrlVar.ThickMin) ;
                    testCase.assertNotEmpty(I,"the synthetic state has no thin ice to constrain")
                    S.BCs.hPosNode=I ; S.BCs.hPosValue=S.CtrlVar.ThickMin+zeros(numel(I),1) ;
                case "ResetThicknessToMinThickness"
                    S.CtrlVar.ResetThicknessToMinThickness=1 ;
                    I=S.F1.h<2*S.CtrlVar.ThickMin ;
                    testCase.assertTrue(any(I),"the synthetic state has no thin ice to reset")
                    S.F1.h(I)=0.7*S.CtrlVar.ThickMin ;   % some thicknesses below ThickMin, so that the reset applies
                case "Enforce_bAboveB"
                    S.CtrlVar.Enforce_bAboveB=true ;
                case "MassBalanceGeometryFeedback"
                    % uvhMatrixAssembly then calls DefineMassBalance.m (via GetMassBalance.m), from the working folder: a temporary
                    % working folder with a DefineMassBalance.m whose mass balance depends on the thickness
                    testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture) ;
                    writelines(["function [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,F)" ; ...
                        "as=0.3+0.001*F.h ; dasdh=0.001+0*F.h ;" ; "ab=-0.002*F.h ; dabdh=-0.002+0*F.h ;" ; "end"],"DefineMassBalance.m") ;
                    rehash ;
                    S.CtrlVar.MassBalanceGeometryFeedback=3 ;
            end
            testCase.verifySameAssembly(S,"option "+Option)
        end

        function testElementTypes(testCase,TriNodes)
            if TriNodes==3
                S=testCase.Base ;
            else
                S=SyntheticUvhState(TriNodes,40,[]) ;
            end
            testCase.verifySameAssembly(S,TriNodes+"-node elements, "+S.MUA.nip+" integration points")
        end

        function testQuadratureRules(testCase,QuadratureDegree)
            S=SyntheticUvhState(3,20,QuadratureDegree) ;
            testCase.verifySameAssembly(S,"linear elements, quadrature degree "+QuadratureDegree+" ("+S.MUA.nip+" points)")
        end

        function testBasalDragAndViscosityOnly(testCase)
            S=testCase.Base ;
            C=S.CtrlVar ; C.OnlyCalcBasalDragAndEffectiveViscosity=true ; C.MassBalanceGeometryFeedback=0 ;
            C.uvhMatrixAssembly.CompiledKernel=true ;
            testCase.assertCompiledPathUsed(C,S.MUA)
            out=cell(2,8) ;
            for k=1:2
                C.uvhMatrixAssembly.CompiledKernel=(k==2) ;
                [out{k,:}]=uvhMatrixAssembly(S.UserVar,S.RunInfo,C,S.MUA,S.F1,S.F1,UaLagrangeVariables,BoundaryConditions()) ;
            end
            testCase.verifyEmpty(out{2,3},"R should be empty") ; testCase.verifyEmpty(out{2,4},"K should be empty")
            Names=["taux","tauy","eta","He"] ;
            for j=5:8
                testCase.verifyLessThanOrEqual(RelDiff(out{1,j},out{2,j}),testCase.Tolerance,Names(j-4)+" at the integration points")
            end
        end

        function testFallbackOutsideCompiledRange(testCase)
            S=SyntheticUvhState(3,20,14) ;   % 42 integration points, beyond the compiled range of 37
            C=S.CtrlVar ; C.uvhMatrixAssembly.CompiledKernel=true ;
            [txt,ok]=evalc("CompiledAssemblyCheck(C,S.MUA,""uvh"")") ;
            testCase.verifyFalse(ok,"the compiled path should not be used for "+S.MUA.nip+" integration points")
            testCase.verifyTrue(contains(txt,"outside the compiled range") || txt=="", ...   % "" if the message was printed earlier in this session
                "unexpected message: "+txt)
            Ca=C ; Ca.uvhMatrixAssembly.CompiledKernel=false ;
            [~,~,Ra,Ka]=uvhMatrixAssembly(S.UserVar,S.RunInfo,Ca,S.MUA,S.F0,S.F1,S.l,S.BCs) ;
            evalc("[~,~,Rb,Kb]=uvhMatrixAssembly(S.UserVar,S.RunInfo,C,S.MUA,S.F0,S.F1,S.l,S.BCs);") ;
            testCase.verifyEqual(Rb,Ra,"after the fallback, R should be identical to the vectorised assembly")
            testCase.verifyEqual(Kb,Ka,"after the fallback, K should be identical to the vectorised assembly")
        end

        function testNewtonSolve(testCase)
            % one implicit uvh time step (Newton-Raphson, line search, active set) with both assemblies, on a well-posed problem:
            % velocities fixed on the upstream and the lateral boundaries
            S=SyntheticUvhState(3,20,[]) ;
            S.CtrlVar.dt=0.001 ;
            x=S.MUA.coordinates(:,1) ; y=S.MUA.coordinates(:,2) ;
            bn=find(x<1 | abs(abs(y)-25e3)<1) ;
            S.BCs.ubFixedNode=bn ; S.BCs.ubFixedValue=S.F1.ub(bn) ; S.BCs.vbFixedNode=bn ; S.BCs.vbFixedValue=S.F1.vb(bn) ;
            C=S.CtrlVar ; C.uvhMatrixAssembly.CompiledKernel=true ;
            testCase.assertCompiledPathUsed(C,S.MUA)
            F=cell(1,2) ; RI=cell(1,2) ;
            for k=1:2
                C.uvhMatrixAssembly.CompiledKernel=(k==2) ;
                evalc("[~,RI{k},F{k}]=uvh(S.UserVar,UaRunInfo,C,S.MUA,S.F0,S.F1,S.l,S.l,S.BCs);") ;
            end
            testCase.verifyTrue(all(RI{2}.Forward.uvhConverged(1)==1),"the uvh solve with the compiled assembly did not converge")
            testCase.verifyEqual(RI{2}.Forward.uvhIterations(1),RI{1}.Forward.uvhIterations(1),"number of Newton iterations")
            for q=["h","ub","vb"]
                testCase.verifyLessThanOrEqual(RelDiff(F{1}.(q),F{2}.(q)),testCase.ToleranceSolve,q+" after one time step")
            end
        end

    end

    methods (Access=private)

        function assertCompiledPathUsed(testCase,CtrlVar,MUA)
            CtrlVar.uvhMatrixAssembly.CompiledKernel=true ;
            [txt,ok]=evalc("CompiledAssemblyCheck(CtrlVar,MUA,""uvh"")") ;
            testCase.assertTrue(ok,"the compiled path is not used, so the comparison would test nothing. "+strtrim(string(txt)))
        end

        function verifySameAssembly(testCase,S,Label)
            % compiled against vectorised assembly, full and residual only
            testCase.assertCompiledPathUsed(S.CtrlVar,S.MUA)
            n=S.MUA.Nnodes ;
            for Ronly=[false true]
                C=S.CtrlVar ; C.uvhMatrixAssembly.Ronly=Ronly ;
                C.uvhMatrixAssembly.CompiledKernel=false ;
                [~,~,Ra,Ka]=uvhMatrixAssembly(S.UserVar,S.RunInfo,C,S.MUA,S.F0,S.F1,S.l,S.BCs) ;
                C.uvhMatrixAssembly.CompiledKernel=true ;
                [~,~,Rb,Kb]=uvhMatrixAssembly(S.UserVar,S.RunInfo,C,S.MUA,S.F0,S.F1,S.l,S.BCs) ;
                What=Label+", "+ternary(Ronly,"residual only","residual and Jacobian")+": " ;
                testCase.verifyFalse(any(isnan(Rb)),What+"NaN in R")
                testCase.verifyLessThanOrEqual(RelDiff(Ra(1:2*n),Rb(1:2*n)),testCase.Tolerance,What+"R, momentum rows")
                testCase.verifyLessThanOrEqual(RelDiff(Ra(2*n+1:end),Rb(2*n+1:end)),testCase.Tolerance,What+"R, thickness rows")
                if ~Ronly
                    testCase.verifyLessThanOrEqual(RelDiff(Ka,Kb),testCase.Tolerance,What+"K")
                end
            end
        end

    end

end

%% local functions

function d=RelDiff(a,b)
% largest difference, relative to the largest absolute value of a
d=full(max(abs(a(:)-b(:))))/max(full(max(abs(a(:)))),realmin) ;
end

function s=ternary(c,a,b)
if c ; s=a ; else ; s=b ; end
end

function S=SyntheticUvhState(TriNodes,nx,QuadratureDegree)
% Small synthetic state for testing the uvh assembly, on a 100 km x 50 km rectangle with nx x nx/2 cells (two triangles each):
% bed from +400 m to -800 m (grounded upstream, floating downstream), thickness from 1000 m to 200 m with a patch thinning to about
% ThickMin, a level-set field crossing the domain near x=90 km, melange fields, and velocities that are not a solution. F0 differs
% slightly from F1 (thickness, velocities, mass balance), as at a time step.
CtrlVar=Ua2D_DefaultParameters() ;
CtrlVar.TriNodes=TriNodes ;
CtrlVar.QuadratureRuleDegree=QuadratureDegree ;
CtrlVar.dt=0.01 ; CtrlVar.time=0 ; CtrlVar.CurrentRunStepNumber=1 ;
CtrlVar.Parallel.uvhAssembly.spmd.isOn=false ; CtrlVar.Parallel.uvAssembly.spmd.isOn=false ;
CtrlVar.uvhMatrixAssembly.ZeroFields=false ; CtrlVar.uvhMatrixAssembly.Ronly=false ; CtrlVar.uvhMatrixAssembly.CompiledKernel=false ;
CtrlVar.InfoLevel=0 ; CtrlVar.InfoLevelNonLinIt=0 ; CtrlVar.doplots=0 ;
CtrlVar.MUA.DecomposeMassMatrix=false ;   % not needed for the uvh assembly (and not possible with a 1-point quadrature rule)

L=100e3 ; W=50e3 ; ny=nx/2 ;
[X,Y]=ndgrid(linspace(0,L,nx+1),linspace(-W/2,W/2,ny+1)) ;
coordinates=[X(:) Y(:)] ;
id=reshape(1:numel(X),nx+1,ny+1) ;
n1=id(1:nx,1:ny) ; n2=id(2:nx+1,1:ny) ; n3=id(2:nx+1,2:ny+1) ; n4=id(1:nx,2:ny+1) ;
connectivity=[n1(:) n3(:) n2(:) ; n1(:) n4(:) n3(:)] ;   % node order as Ua expects
if TriNodes~=3
    evalc("[coordinates,connectivity]=ChangeElementType(coordinates,connectivity,TriNodes);") ;
end
evalc("MUA=CreateMUA(CtrlVar,connectivity,coordinates);") ;

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ; n=MUA.Nnodes ; xi=x/L ; eta=y/W ;
F1=UaFields ; F1.x=x ; F1.y=y ;
F1.rho=917+zeros(n,1) ; F1.rhow=1030 ; F1.g=9.81/1000 ; F1.alpha=0 ;
F1.S=zeros(n,1) ;
F1.B=400-1200*xi+100*cos(2*pi*eta) ;
h=1000*(1-0.8*xi)+50*sin(3*pi*xi).*cos(2*pi*eta) ;
thin=min(1,2*max(0,1-((x-0.9*L).^2+(y-0.25*W).^2)/(10e3)^2)) ;    % a patch where the ice thins to about ThickMin
F1.h=h.*(1-thin)+thin.*(CtrlVar.ThickMin*(1.3+0.5*sin(5*pi*xi).^2)) ;
[F1.b,F1.s]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow) ;
F1.ub=50+500*xi.^2+20*sin(2*pi*eta) ; F1.vb=30*sin(pi*xi).*cos(pi*eta) ; F1.ud=zeros(n,1) ; F1.vd=zeros(n,1) ;
F1.AGlen=AGlenVersusTemp(-10)*(1+0.3*xi) ; F1.n=3+zeros(n,1) ;
F1.C=0.01*(1+0.5*sin(2*pi*xi)) ; F1.m=3+zeros(n,1) ; F1.q=1+zeros(n,1) ; F1.muk=0.5+zeros(n,1) ; F1.V0=100+zeros(n,1) ;
F1.as=0.3+zeros(n,1) ; F1.ab=-1*(F1.h<F1.rhow*(F1.S-F1.B)./F1.rho) ; F1.dasdh=zeros(n,1) ; F1.dabdh=zeros(n,1) ;
F1.LSF=0.9*L-x+5e3*sin(2*pi*eta) ;
F1.uo=0.3*F1.ub ; F1.vo=0.3*F1.vb ; F1.ua=5*F1.ub+10 ; F1.va=5*F1.vb-10 ;
F1.Co=10*F1.C ; F1.mo=3+zeros(n,1) ; F1.Ca=100*F1.C ; F1.ma=3+zeros(n,1) ;
F0=F1 ; F0.h=F1.h.*(1-0.01*cos(pi*xi)) ; F0.ub=0.98*F1.ub ; F0.vb=0.97*F1.vb ; F0.as=0.9*F1.as ;
evalc("CtrlVar=SetSmoothFloorWidths(CtrlVar,F1);") ;

S=struct("UserVar",[],"RunInfo",UaRunInfo,"CtrlVar",CtrlVar,"MUA",MUA,"F0",F0,"F1",F1,"l",UaLagrangeVariables,"BCs",BoundaryConditions()) ;
end
