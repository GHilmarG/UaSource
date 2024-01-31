





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

        else  % this is the otherwise default sequencial assembly

            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F) ;

        end

        %% Test spmd speedup

        if CtrlVar.Parallel.isTest  && ~CtrlVar.uvMatrixAssembly.Ronly


            % There might be an argument for doing this a few times, but if one does a number of NR iterations, which generally is the
            % case, then that should not be needed.
            tSeq=tic ;  [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F); tSeq=toc(tSeq) ;


            if CtrlVar.Parallel.uvAssembly.spmd.isOn

                tSPMD=tic ;  [RuvSPMD,KuvSPMD,Tint,Fext]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F); tSPMD=toc(tSPMD) ;
                fprintf(' tSeq=%f sec     tSPMD=%f sec \t MUA.Nnodes=%i \t     norm(Ruv-RuvSPMD)/norm(Ruv)=%g \t     norm(diag(Kuv)-diag(KuvSPMD))/norm(diag(Kuv))=%g \n',...
                    tSeq,tSPMD,MUA.Nnodes,norm(full(Ruv-RuvSPMD))/norm(full(Ruv)),norm(diag(Kuv)-diag(KuvSPMD))/norm(diag(Kuv)))
            end

            if CtrlVar.Parallel.uvAssembly.parfeval.isOn
                tParfeval=tic ;  [RuvParfeval,KuvParfeval,Tint,Fext]=uvMatrixAssemblySSTREAM_Parfeval(CtrlVar,MUA,F); tParfeval=toc(tParfeval) ;
                fprintf(' tSeq=%f sec  tParfval=%f sec \t MUA.Nnodes=%i \t norm(Ruv-RuvParfeval)/norm(Ruv)=%g \t norm(diag(Kuv)-diag(KuvParfeval))/norm(diag(Kuv))=%g \n',...
                    tSeq,tParfeval,MUA.Nnodes,norm(full(Ruv-RuvParfeval))/norm(full(Ruv)),norm(diag(Kuv)-diag(KuvParfeval))/norm(diag(Kuv)))
            end

            fprintf(" \n [------------------- Comparing speedup in uv assembly using parallel options: \n ")



            if CtrlVar.Parallel.uvAssembly.spmd.isOn
                fprintf(' Speedup in matrix assembly by using the parallel SPMD approach is %f \n',tSeq/tSPMD)
            end

            if CtrlVar.Parallel.uvAssembly.parfeval.isOn
                fprintf(' Speedup in matrix assembly by using the parallel parfeval approach is %f \n',tSeq/tParfeval)
            end

            fprintf(" -------------------------------------------------------------------------------------] \n ")

        end



    case "sstream-rho"

        [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAMrho(CtrlVar,MUA,F) ;

    otherwise

        error("UA:uvMatrixAssemblyCaseNotFound","Case Not Found")

end


RunInfo.CPU.Assembly.uv=RunInfo.CPU.Assembly.uv+toc(tAssembly);

end



