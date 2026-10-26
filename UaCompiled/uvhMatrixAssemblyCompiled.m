function [R,K,tauxInt,tauyInt,etaInt,HeInt]=uvhMatrixAssemblyCompiled(CtrlVar,MUA,F0,F1,BCs1,LSFMask,dadh,Ronly,RK,KernelFull,KernelRonly)

%% uvh assembly with the compiled element-based kernel
%
% Called from uvhMatrixAssembly.m if CtrlVar.uvhMatrixAssembly.CompiledKernel=true and uvhCompiledKernelCheck.m allows it (RK, KernelFull
% and KernelRonly come from there). Returns the residual R (=Tint-Fext) and, unless Ronly, the Jacobian K, as uvhMatrixAssembly.m does.
% Basal drag, effective viscosity and floating mask at the integration points are only calculated if requested (nargout>2).
%
% The kernel (uvhElementKernel.m) gathers the nodal values itself from the nodal vectors and the connectivity, and sums the element-matrix
% values directly into the nonzeros of K (in CSC order, see uvhCSCPattern.m); K is then created directly from its CSC arrays
% (SparseFromCSC). The residual is summed with the precomputed indices of the assembly pattern (P.iR). See README.md in this folder.
% (25 Oct 2026)
%%

Nn=MUA.Nnodes; neq=3*Nn; P=MUA.uvhAssemblyPattern;
IntOut = nargout>2 ;
col=@(x) reshape(double(x),[],1);    % nodal column vector (an empty field gives a 0x1 vector)
hBCsMask=zeros(Nn,1); hBCsMask(BCs1.hFixedNode)=1; hBCsMask(BCs1.hPosNode)=1;
if CtrlVar.IncludeMelangeModelPhysics
    uo=col(F1.uo); vo=col(F1.vo); Co=col(F1.Co); mo=col(F1.mo); ua=col(F1.ua); va=col(F1.va); Ca=col(F1.Ca); ma=col(F1.ma);
else
    uo=zeros(0,1); vo=uo; Co=uo; mo=uo; ua=uo; va=uo; Ca=uo; ma=uo;
end
if Ronly
    Cp=[]; invPtr=zeros(0,1,'uint32'); invIdx=zeros(0,1,'uint32');
else
    Cp=uvhCSCPattern(P); invPtr=Cp.invPtr; invIdx=Cp.invIdx;
end
% Deriv and DetJ are passed as separate arrays: fields of a structure input would be copied into the MEX file, arrays are not
MUAk=struct("Nele",MUA.Nele,"nod",MUA.nod,"nip",MUA.nip,"points",MUA.points,"weights",MUA.weights(:).',"EleAreas",MUA.EleAreas(:));
args={MUAk,MUA.Deriv,MUA.DetJ,MUA.connectivity,col(F1.b),col(F1.h),col(F1.ub),col(F1.vb),col(F1.AGlen),col(F1.n),col(F1.C),col(F1.m), ...
      col(F1.q),col(F1.muk),col(F1.V0),col(F0.h),col(F0.ub),col(F0.vb),col(F0.as),col(F0.ab),col(F1.as),col(F1.ab),col(dadh),col(F1.B),col(F1.S), ...
      col(F1.rho),col(LSFMask),hBCsMask,uo,vo,Co,mo,ua,va,Ca,ma,RK,F1.rhow,F1.g,cos(F1.alpha),sin(F1.alpha),CtrlVar.dt,IntOut,invPtr,invIdx};
if Ronly
    [Tx,Fx,Ty,Fy,Th,Fh,~,tauxInt,tauyInt,etaInt,HeInt]=feval(KernelRonly,args{:});
else
    [Tx,Fx,Ty,Fy,Th,Fh,vs,tauxInt,tauyInt,etaInt,HeInt]=feval(KernelFull,args{:});
end

K=[];
if isfield(CtrlVar,"OnlyCalcBasalDragAndEffectiveViscosity") && CtrlVar.OnlyCalcBasalDragAndEffectiveViscosity
    R=[]; return
end

% residual, summed with the precomputed indices (same order as in uvhMatrixAssembly.m, so the same result)
Tval=reshape(permute(cat(3,Tx,Ty,Th),[1 3 2]),[],1);
Fval=reshape(permute(cat(3,Fx,Fy,Fh),[1 3 2]),[],1);
Tint=accumarray(P.iR,Tval,[neq 1]);
Fext=accumarray(P.iR,Fval,[neq 1]);
R=Tint-Fext;
if any(isnan(R)) ||  any(isnan(Tint)) || any(isnan(Fext))
    fprintf("nan in R or Tint or Fext")
end

if ~Ronly
    if any(isnan(vs(Cp.diagIdx)))
        error(' NaN in K ' ) ;
    end
    K=SparseFromCSC(Cp.jc,Cp.ir,vs,neq);
end

end
