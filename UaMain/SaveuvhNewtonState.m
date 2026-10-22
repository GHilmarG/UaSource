function SaveuvhNewtonState(UserVar,CtrlVar,MUA,F0,F1,l1,BCs1,duvh,dl,iteration,r0)

%% Saves the state of a uvh Newton-Raphson iteration, for testing the Jacobian afterwards (TestuvhJacobian.m)
%
%   SaveuvhNewtonState(UserVar,CtrlVar,MUA,F0,F1,l1,BCs1,duvh,dl,iteration,r0)
%
% (15 Oct 2026) Called in each uvh Newton-Raphson iteration if CtrlVar.uvhSaveNewtonSystem.Save=true. The state is saved if the
% iteration is selected (see uvhNewtonStateSelected.m), and if fewer than CtrlVar.uvhSaveNewtonSystem.MaxStates states have been saved.
%
% All states of a run are saved to one file, CtrlVar.uvhSaveNewtonSystem.FileName, each as a separate variable State_0001,
% State_0002, ... The state contains the inputs of uvhAssembly (UserVar, CtrlVar, F0, F1, l1, BCs1), from which the residual and the
% Jacobian are recalculated in the test, the Newton direction of the iteration (duvh, dl), and the context (run step, time, time step,
% solve number within the run step, Newton-Raphson iteration, r0). The mesh (MUA) is saved only once, as MUA_001, and again only if it
% changes (MUA_002, ...); each state refers to its mesh (MeshIndex). If the file exists already, new states are appended.
%
%%

persistent FileName nSaved MeshKey MeshIndex LastStep SolveIndex

Opt=CtrlVar.uvhSaveNewtonSystem ;
fn=string(Opt.FileName) ;

%% new file (or first call): initialise, continuing the numbering of an existing file
if isempty(FileName) || fn~=FileName
    FileName=fn ; nSaved=0 ; MeshKey=[] ; MeshIndex=0 ; LastStep=NaN ; SolveIndex=0 ;
    if isfile(fn)
        w=whos("-file",fn) ; names=string({w.name}) ;
        nSaved=sum(startsWith(names,"State_")) ;
        MeshIndex=sum(startsWith(names,"MUA_")) ;
    end
end

%% count the solves within a run step (a solve starts with Newton iteration 1)
if CtrlVar.CurrentRunStepNumber~=LastStep
    LastStep=CtrlVar.CurrentRunStepNumber ; SolveIndex=0 ;
end
if iteration==1
    SolveIndex=SolveIndex+1 ;
end

%% is this iteration to be saved?
MaxStates=100 ; if isfield(Opt,"MaxStates") ; MaxStates=Opt.MaxStates ; end
if ~uvhNewtonStateSelected(Opt,CtrlVar,iteration) || nSaved>=MaxStates
    return
end

%% the mesh, if new or changed
Key=[MUA.Nnodes MUA.Nele sum(MUA.coordinates(:)) sum(MUA.connectivity(:))] ;
NewMesh = ~isequal(Key,MeshKey) ;
Out=struct() ;
if NewMesh
    MeshKey=Key ; MeshIndex=MeshIndex+1 ;
    MUAsave=MUA ;
    if isfield(MUAsave,"dM") ; MUAsave.dM=[] ; end            % decompositions can not be saved
    if isfield(MUAsave,"workers") ; MUAsave.workers=[] ; end
    Out.(sprintf("MUA_%03d",MeshIndex))=MUAsave ;
end

%% the state
nSaved=nSaved+1 ;
S.UserVar=UserVar ; S.CtrlVar=CtrlVar ; S.F0=F0 ; S.F1=F1 ; S.l1=l1 ; S.BCs1=BCs1 ;
S.duvh=duvh ; S.dl=dl ;
S.RunStep=CtrlVar.CurrentRunStepNumber ; S.time=CtrlVar.time ; S.dt=CtrlVar.dt ;
S.Solve=SolveIndex ; S.Iteration=iteration ; S.r0=r0 ; S.MeshIndex=MeshIndex ;
Out.(sprintf("State_%04d",nSaved))=S ;

if isfile(fn)
    save(fn,"-struct","Out","-append") ;
else
    save(fn,"-struct","Out","-v7.3") ;
end

if isfield(CtrlVar,"InfoLevelNonLinIt") && CtrlVar.InfoLevelNonLinIt>=1
    fprintf(" uvhSaveNewtonSystem: state %i saved to %s (run step %i, solve %i, Newton iteration %i, r0=%g)\n",...
        nSaved,fn,S.RunStep,S.Solve,S.Iteration,r0) ;
end

end