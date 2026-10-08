function [UserVar,RunInfo,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)
    
    narginchk(8,8)

    nargoutchk(3,4)
    tAssembly=tic;
    
    if CtrlVar.Parallel.uvhAssembly.spmd.isOn
        [UserVar,RunInfo,R,K]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
    else
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
    end
    %% Modification (8 Oct 2026): lumped nodal thickness penalty ("Variant 2").
    % The penalty is added as a nodal term, using lumped masses m_i (row sums of the mass matrix), and depends only on
    % the nodal thickness h_i. Nodes with h boundary conditions or active thickness constraints are excluded (w_i=0).
    % Residual and Jacobian contributions are diagonal. Not added to the ZeroFields (normalisation) reference.
    % Activated by CtrlVar.ThicknessPenalty=true together with CtrlVar.ThicknessPenaltyNodal=true, in which case the
    % integration-point penalty in uvhAssemblyIntPointImplicitSUPG_v2 is switched off.
    if isfield(CtrlVar,"ThicknessPenalty") && CtrlVar.ThicknessPenalty ...
            && isfield(CtrlVar,"ThicknessPenaltyNodal") && CtrlVar.ThicknessPenaltyNodal ...
            && ~CtrlVar.uvhMatrixAssembly.ZeroFields
        Nn=MUA.Nnodes;
        if isfield(MUA,"M") && ~isempty(MUA.M) ; mLump=MUA.M*ones(Nn,1) ; else ; mLump=MassMatrix2D1dof(MUA)*ones(Nn,1) ; end
        wPen=ones(Nn,1) ; wPen(BCs1.hFixedNode)=0 ; wPen(BCs1.hPosNode)=0 ;
        [aPen,daPendh]=ThicknessPenaltyMassBalanceFeedback(CtrlVar,F1.h) ;
        cPen=CtrlVar.dt*CtrlVar.theta*F1.rho(:).*mLump.*wPen ;
        ih=(2*Nn+1:3*Nn)' ;
        R(ih)=R(ih)-cPen.*aPen ;
        if ~CtrlVar.uvhMatrixAssembly.Ronly && nargout>=4 && ~isempty(K)
            K=K+sparse(ih,ih,-cPen.*daPendh,size(K,1),size(K,2)) ;
        end
    end

    tAssembly=toc(tAssembly);
    RunInfo.CPU.Assembly.uvh=RunInfo.CPU.Assembly.uvh+tAssembly;

    if CtrlVar.InfoLevelCPU>=10
        if CtrlVar.Parallel.uvhAssembly.spmd.isOn
            if CtrlVar.uvhMatrixAssembly.Ronly
                fprintf("SPMD Force vector assembly in %g sec \n",tAssembly)
            else
                fprintf("SPMD Matrix assembly in %g sec \n",tAssembly)
            end
        else
            if CtrlVar.uvhMatrixAssembly.Ronly
                fprintf("Force vector assembly in %g sec \n",tAssembly)
            else
                fprintf("Matrix assembly in %g sec \n",tAssembly)
            end
        end
    end



    %%

    if CtrlVar.Parallel.isTest && ~CtrlVar.uvhMatrixAssembly.Ronly

        tSeq=tic;
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
        tSeq=toc(tSeq) ;

        CtrlVar.Parallel.uvhAssembly.spmd.isOn=true ;
        MUA=UpdateMUA(CtrlVar,MUA);  % just in case workers have not been set up
        tSPMD=tic;
        [UserVar,RunInfo,Rspmd,Kspmd]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
        tSPMD=toc(tSPMD);


        if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=0;
            else
                CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;
            end
        end

        fprintf('\n ----------------------------- Info on parallel uvh SPMD assembly performance : Nele=%i  \t nWorkers=%i \n',MUA.Nele,CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
        fprintf('#Ele=%i \t SPMD used for uvh assembly:  tSeq=%f \t tSPMD=%f \t speedup=%g \n',MUA.Nele,tSeq,tSPMD,tSeq/tSPMD) ;
        fprintf(' R-Rspmd=%g \t K-Kspmd=%g   \n',full(norm(R-Rspmd)/norm(R)),normest(K-Kspmd)/normest(K))
        fprintf(' ----------------------------- \n')




    end
    %

end
