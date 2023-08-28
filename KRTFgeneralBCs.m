function [Ruv,Kuv,Tint,Fext]=KRTFgeneralBCs(CtrlVar,MUA,F,ZeroFields)

%%
%
%   [Ruv,Kuv,Tint,Fext]=KRTFgeneralBCs(CtrlVar,MUA,F,ZeroFields)
%
%  uv SSA/SSTREAM assembly
%
%
%
%%

if nargin<4
    CtrlVar.uvAssembly.ZeroFields=false;
else
    CtrlVar.uvAssembly.ZeroFields=ZeroFields;
end


if nargout==1
    
    CtrlVar.uvMatrixAssembly.Ronly=1;
else
    CtrlVar.uvMatrixAssembly.Ronly=0;
end



%[Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(CtrlVar,MUA,F,ZeroFields);

%tAssembly=tic;

if CtrlVar.Parallel.uvAssembly.spmd.isOn && ~CtrlVar.uvMatrixAssembly.Ronly
    [Ruv,Kuv,Tint,Fext]=uvAssemblySPMD(CtrlVar,MUA,F);
else
    [Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(CtrlVar,MUA,F);
end


%RunInfo.CPU.Assembly.uv=toc(tAssembly)+RunInfo.CPU.Assembly.uv;

%%

 
if CtrlVar.Parallel.isTest && CtrlVar.Parallel.uvAssembly.spmd.isOn  && ~CtrlVar.uvMatrixAssembly.Ronly
    
    tSeq=tic;
    [Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(CtrlVar,MUA,F);
    tSeq=toc(tSeq) ;
    
    tSPMD=tic;
    [Ruv2,Kuv,Tint,Fext]=uvAssemblySPMD(CtrlVar,MUA,F);
    tSPMD=toc(tSPMD);
    
    if ~CtrlVar.uvMatrixAssembly.Ronly
        
        fprintf('\n \n ----------------------------- \n')
        fprintf(' tSeq=%f \t tSPMD=%f \t MUA.Nele=%i \t norm(Ruv-Ruv2)/norm(Ruv)=%g \n',tSeq,tSPMD,MUA.Nele,norm(full(Ruv-Ruv2))/norm(full(Ruv)))
        %fprintf(RunInfo.File.fid,' tSeq=%f \t tSPMD=%f \t MUA.Nele=%i \t %f \n',tSeq,tSPMD,MUA.Nele,norm(full(Ruv-Ruv2))/norm(full(Ruv)));
    end
    
end
%

end





