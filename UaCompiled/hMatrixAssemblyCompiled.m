function [f0,K]=hMatrixAssemblyCompiled(CtrlVar,MUA,F0,F1,a0,a1,dadh,LSFMask,hBCsMask,Ronly,RK,KernelFull,KernelRonly)

%% h assembly (mass conservation, velocities given) with the compiled element-based kernel
%
% Called from MassContinuityEquationAssembly.m if CtrlVar.hMatrixAssembly.CompiledKernel=true and the check of the compiled assembly
% allows it (RK, KernelFull and KernelRonly come from there), after the common preparation there:
%
%   a0, a1   : mass balance (as+ab) at the beginning and end of the time step (a1 updated for the mass-balance geometry feedback)
%   dadh     : nodal derivative of the mass balance with respect to h (zero without mass-balance geometry feedback)
%   LSFMask  : level-set mask F1.LSFMask.NodesOut (only used with the level-set mass-balance feedback; may be empty otherwise)
%   hBCsMask : nodes with thickness boundary conditions or active thickness constraints (logical or 0/1, Nnodes x 1)
%
% Returns the residual f0 (full vector) and, unless Ronly, the matrix K (K dh = -f0), as the vectorised assembly in
% MassContinuityEquationAssembly.m does, before the lumped nodal thickness penalty, which is added there for both.
%
% The kernel (hElementKernel.m) gathers the nodal values itself from the nodal vectors and the connectivity, sums the element residuals
% into f0 and the element-matrix values into the nonzeros of K (in CSC order, see CompiledAssemblyPattern.m); K is then created directly
% from its CSC arrays (SparseFromCSC). Complex-valued inputs are errors (the kernel is compiled for real values).
%
% (31 Oct 2026)
%%

Nn=MUA.Nnodes; neq=Nn; P=MUA.hAssemblyPattern;
col=@(x) reshape(double(x),[],1);    % nodal column vector (an empty input gives a 0x1 vector)

Cp=CompiledAssemblyPattern(P);   % per mesh, cached: for the summation of the residual (and, unless Ronly, of the matrix) in the kernel

Fields={F0.h,F1.h,F0.ub,F0.vb,F1.ub,F1.vb,F1.rho,a0,a1,dadh,LSFMask,hBCsMask} ;
if ~all(cellfun(@isreal,Fields))
    error("hMatrixAssemblyCompiled:ComplexInput","Complex-valued input to the h assembly.")
end

% Deriv and DetJ are passed as separate arrays: fields of a structure input would be copied into the MEX file, arrays are not
MUAk=struct("Nele",MUA.Nele,"nod",MUA.nod,"nip",MUA.nip,"points",MUA.points,"weights",MUA.weights(:).',"EleAreas",MUA.EleAreas(:));
args=[{MUAk,MUA.Deriv,MUA.DetJ,MUA.connectivity}, cellfun(col,Fields,'UniformOutput',false), {RK,Cp.rPtr,Cp.rIdx,Cp.invPtr,Cp.invIdx}];

% A failure of the MEX file itself (eg built for another MATLAB release or another platform, or with an outdated interface) is reported
% as UaCompiled:KernelFailed: the vectorised assembly is then used instead.
try
    if Ronly
        f0=feval(KernelRonly,args{:});
    else
        [f0,vs]=feval(KernelFull,args{:});
    end
catch ME
    if Ronly ; Kernel=KernelRonly ; else ; Kernel=KernelFull ; end
    E=MException("UaCompiled:KernelFailed","%s: %s",Kernel,ME.message) ;
    throw(addCause(E,ME)) ;
end

K=[];
if ~Ronly
    try
        K=SparseFromCSC(Cp.jc,Cp.ir,vs,neq);
    catch ME
        E=MException("UaCompiled:KernelFailed","SparseFromCSC: %s",ME.message) ;
        throw(addCause(E,ME)) ;
    end
end

end
