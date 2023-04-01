function [Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(CtrlVar,MUA,F)

%
% Ruv=Tint-Fext;
% Tint   : internal nodal forces
% Fint   : external nodal forces

narginchk(3,5)
nargoutchk(1,5)

switch lower(CtrlVar.FlowApproximation)

    case "sstream"


        if CtrlVar.Parallel.uvAssembly.spmdInt.isOn  && ~CtrlVar.uvMatrixAssembly.Ronly

            % This uses spmd over integration points with one worker per integration point
            % Does speed things up somewhat with increasing number of elements and nip
            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAMspmd(CtrlVar,MUA,F) ;

        else

            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F) ;

        end

        %% Do both sequential and parallel spmd assembly and compare times, and save data for further analysis
        if CtrlVar.Parallel.isTest && CtrlVar.Parallel.uvAssembly.spmdInt.isOn &&  ~CtrlVar.uvMatrixAssembly.Ronly

            tSeq=tic;
            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F) ;
            tSeq=toc(tSeq) ;

            tSPMD=tic;
            [Ruv2,Kuv2,Tint2,Fext2]=uvMatrixAssemblySSTREAMspmd(CtrlVar,MUA,F) ;
            tSPMD=toc(tSPMD);

            fprintf('\n \n ----------------------------- \n')
            fprintf(' tSeq=%f \t tSPMD=%f \t speedup=%f \n',tSeq,tSPMD,tSeq/tSPMD)
            fprintf(' norm(Ruv-Ruv2)/norm(Ruv)=%g \t Kuv error=%g \n',norm(full(Ruv-Ruv2))/norm(full(Ruv)),norm(diag(Kuv2-Kuv))/max(abs(diag(Kuv))))
             
            

            FileName="TestSaveuvMatrixAssembly"+"nEle"+num2str(MUA.Nele ...
                ) ;
            fprintf("Saving MUA and F in %s \n",FileName)
            save(FileName,"CtrlVar","MUA","F")

        end



    case "sstream-rho"

        [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAMrho(CtrlVar,MUA,F) ;

    otherwise

        error("UA:uvMatrixAssemblyCaseNotFound","Case Not Found")

end


end



