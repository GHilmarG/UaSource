function [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblyCompiled(CtrlVar,MUA,F,Ronly,RK,KernelFull,KernelRonly)

%% uv assembly with the compiled element-based kernel
%
% Called from uvMatrixAssemblySSTREAM.m if CtrlVar.uvMatrixAssembly.CompiledKernel=true and the check of the compiled assembly allows it
% (RK, KernelFull and KernelRonly come from there), after the common preparation (ZeroFields, checks). Returns, as
% uvMatrixAssemblySSTREAM.m does, the residual Ruv=Tint-Fext (full vectors) and, unless Ronly, the symmetrised tangent matrix Kuv.
%
% The kernel (uvElementKernel.m) gathers the nodal values itself from the nodal vectors and the connectivity, sums the element-matrix
% values directly into the nonzeros of Kuv (in CSC order, see CompiledAssemblyPattern.m) and symmetrises them; Kuv is then created
% directly from its CSC arrays (SparseFromCSC). The residual is summed with the precomputed indices of the assembly pattern (P.iR), in the kernel (since 31 Oct 2026), in
% the same order as in uvMatrixAssemblySSTREAM.m. See README.md in this folder.
%
% Complex-valued inputs are errors (the kernel is compiled for real values).
%
% (30 Oct 2026)
%%

Nn=MUA.Nnodes; neq=2*Nn; P=MUA.uvAssemblyPattern;
col=@(x) reshape(double(x),[],1);    % nodal column vector (an empty field gives a 0x1 vector)

if CtrlVar.IncludeMelangeModelPhysics
    uo=col(F.uo); vo=col(F.vo); Co=col(F.Co); mo=col(F.mo); ua=col(F.ua); va=col(F.va); Ca=col(F.Ca); ma=col(F.ma);
else
    uo=zeros(0,1); vo=uo; Co=uo; mo=uo; ua=uo; va=uo; Ca=uo; ma=uo;
end
Cp=CompiledAssemblyPattern(P);   % per mesh, cached: for the summation of the residual (and, unless Ronly, of the matrix) in the kernel

Fields={F.s,F.h,F.ub,F.vb,F.AGlen,F.n,F.C,F.m,F.q,F.muk,F.V0,F.B,F.S,F.rho,uo,vo,Co,mo,ua,va,Ca,ma} ;
if ~all(cellfun(@isreal,Fields)) || ~isreal(F.rhow) || ~isreal(F.g) || ~isreal(F.alpha)
    error("uvMatrixAssemblyCompiled:ComplexInput","Complex-valued input to the uv assembly.")
end

% Deriv and DetJ are passed as separate arrays: fields of a structure input would be copied into the MEX file, arrays are not
MUAk=struct("Nele",MUA.Nele,"nod",MUA.nod,"nip",MUA.nip,"points",MUA.points,"weights",MUA.weights(:).');
args=[{MUAk,MUA.Deriv,MUA.DetJ,MUA.connectivity}, cellfun(col,Fields(1:14),'UniformOutput',false), {uo,vo,Co,mo,ua,va,Ca,ma}, ...
      {RK,F.rhow,F.g,cos(F.alpha),sin(F.alpha),Cp.rPtr,Cp.rIdx,Cp.invPtr,Cp.invIdx,Cp.perm}];

% A failure of the MEX file itself (eg built for another MATLAB release or another platform, or with an outdated interface) is reported
% as UaCompiled:KernelFailed: the vectorised assembly is then used instead.
try
    if Ronly
        [Tint,Fext]=feval(KernelRonly,args{:});
    else
        [Tint,Fext,vs]=feval(KernelFull,args{:});
    end
catch ME
    if Ronly ; Kernel=KernelRonly ; else ; Kernel=KernelFull ; end
    E=MException("UaCompiled:KernelFailed","%s: %s",Kernel,ME.message) ;
    throw(addCause(E,ME)) ;
end

% (31 Oct 2026) Tint and Fext are summed in the kernel, in the same order as accumarray(P.iR,...) in uvMatrixAssemblySSTREAM.m
Ruv=Tint-Fext;

Kuv=[];
if ~Ronly
    try
        Kuv=SparseFromCSC(Cp.jc,Cp.ir,vs,neq);
    catch ME
        E=MException("UaCompiled:KernelFailed","SparseFromCSC: %s",ME.message) ;
        throw(addCause(E,ME)) ;
    end
end

end
