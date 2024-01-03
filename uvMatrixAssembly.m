function [Ruv,Kuv,Tint,Fext,MUAworkers]=uvMatrixAssembly(CtrlVar,MUA,F,MUAworkers)

%
% Ruv=Tint-Fext;
% Tint   : internal nodal forces
% Fint   : external nodal forces

narginchk(4,4)
nargoutchk(1,5)



switch lower(CtrlVar.FlowApproximation)

    case "sstream"

        if CtrlVar.uvMatrixAssembly.Ronly  % this is an assembly of the right-hand-side only, this is fast anyhow

            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F) ;

        else  % This involves assembly of the matrix, this might be worth trying to speed up using parallel options

            if CtrlVar.Parallel.uvAssembly.parfeval.isOn

                [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_Parfeval(CtrlVar,MUA,F);

            elseif CtrlVar.Parallel.uvAssembly.spmd.isOn

              
                [Ruv,Kuv,Tint,Fext,MUAworkers]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,MUAworkers);

            else

                [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F) ;

            end

            if CtrlVar.Parallel.isTest

                % MUAworkers=[];

                tSeq=tic ;  [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F); tSeq=toc(tSeq) ;
                % tSeq2=tic ;  [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F); tSeq2=toc(tSeq2) ;   [tSeq tSeq2]

                if CtrlVar.Parallel.uvAssembly.spmd.isOn
                    MUAworkers=[]; 
                    tSPMD=tic ;  [RuvSPMD,KuvSPMD,Tint,Fext,MUAworkers]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,MUAworkers); tSPMD=toc(tSPMD) ;
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

        end

    case "sstream-rho"

        [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAMrho(CtrlVar,MUA,F) ;

    otherwise

        error("UA:uvMatrixAssemblyCaseNotFound","Case Not Found")

end


end



