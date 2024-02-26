





function [RunInfo,Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F)

%
% Ruv=Tint-Fext;
% Tint   : internal nodal forces
% Fint   : external nodal forces

narginchk(4,4)
nargoutchk(2,4)



if nargout==2 && ~CtrlVar.uvMatrixAssembly.Ronly
    error("KRTFgeneralBCs: More than one output, but assembly only required for R. ")
end

tAssembly=tic;

switch lower(CtrlVar.FlowApproximation)

    case "sstream"


        if CtrlVar.Parallel.uvAssembly.parfeval.isOn

            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_Parfeval(CtrlVar,MUA,F);

        elseif CtrlVar.Parallel.uvAssembly.spmd.isOn


            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F);

        else  % this is the otherwise default sequential assembly

            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F) ;

        end

        %% Test spmd speedup

        if CtrlVar.Parallel.isTest  && CtrlVar.Parallel.uvAssembly.spmd.isOn && ~CtrlVar.uvMatrixAssembly.Ronly


            % There might be an argument for doing this a few times, but if one does a number of NR iterations, which generally is the
            % case, then that should not be needed.
            tSeq=tic ;  [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F); tSeq=toc(tSeq) ;


            if CtrlVar.Parallel.uvAssembly.spmd.isOn

                tSPMD=tic ;  [RuvSPMD,KuvSPMD,Tint,Fext]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F); tSPMD=toc(tSPMD) ;
               
            end

            if CtrlVar.Parallel.uvAssembly.parfeval.isOn
                tParfeval=tic ;  [RuvParfeval,KuvParfeval,Tint,Fext]=uvMatrixAssemblySSTREAM_Parfeval(CtrlVar,MUA,F); tParfeval=toc(tParfeval) ;
                fprintf(' tSeq=%f sec  tParfval=%f sec \t MUA.Nnodes=%i \t norm(Ruv-RuvParfeval)/norm(Ruv)=%g \t norm(diag(Kuv)-diag(KuvParfeval))/norm(diag(Kuv))=%g \n',...
                    tSeq,tParfeval,MUA.Nnodes,norm(full(Ruv-RuvParfeval))/norm(full(Ruv)),norm(diag(Kuv)-diag(KuvParfeval))/norm(diag(Kuv)))
            end

            if CtrlVar.Parallel.uvAssembly.spmd.isOn
                fprintf('\n ----------------------------- Info on parallel uv SPMD assembly performance : nEle=%i  \t nWorkers=%i \n',MUA.Nele,CtrlVar.Parallel.uvAssembly.spmd.nWorkers)
                fprintf('SPMD used for uv assembly:  tSeq=%f \t tSPMD=%f \t Speedup=%g \n',tSeq,tSPMD,tSeq/tSPMD) ;
                fprintf(' R-Rspmd=%g \t K-Kspmd=%g   \n',full(norm(Ruv-RuvSPMD)/norm(Ruv)),normest(Kuv-KuvSPMD)/normest(Kuv))
                fprintf(' ----------------------------- \n')
            end

        end


    case "sstream-rho"

        [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAMrho(CtrlVar,MUA,F) ;

    otherwise

        error("UA:uvMatrixAssemblyCaseNotFound","Case Not Found")

end


RunInfo.CPU.Assembly.uv=RunInfo.CPU.Assembly.uv+toc(tAssembly);

end



