




function [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l)
         

%%
%
% solves for velocities
%
% If we write the finite-element formulation of the forward model on the form
% 
% $$\mathbf{F}_{\mathbf{v}} = \mathbf{0} $$
%
% and then write
%
% $$ \mathbf{f}_{\mathbf{v}} - \mathbf{v} = 0$$
%
%
% then the vector uv is $\mathbf{f}$
%
% the vector Ruv is $\mathbf{F}_{\mathbf{v}}$
%
% and the matrix Kuv = $d_{\mathbf{v}} \mathbf{F}_{\mathbf{uv}}$ 
%
%%

nargoutchk(4,8);
narginchk(7,7)

tdiagnostic=tic;

if ( CtrlVar.Parallel.uvAssembly.spmd.isOn || CtrlVar.Parallel.uvhAssembly.spmd.isOn  )

    poolobj = gcp('nocreate');

    if isempty(poolobj)

        fprintf("uv SPMD assembly is set to true, but parallel pool is empty. \n")
        fprintf(" Create a parallel pool ahead of the call to %ca.\n",218)
    else
        CtrlVar.Parallel.uvAssembly.spmd.nWorkers=poolobj.NumWorkers;
    end

end

if isscalar(F.m)
    F.m=zeros(MUA.Nnodes,1)+F.m;

end

if isscalar(F.n)
    F.n=zeros(MUA.Nnodes,1)+F.n;
end

if isempty(l)
    l=UaLagrangeVariables; 
end



[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);



if CtrlVar.LevelSetMethod &&  ~isnan(CtrlVar.LevelSetDownstreamAGlen)  &&  ~isnan(CtrlVar.LevelSetDownstream_nGlen)
    if isempty(F.LSFMask)  % This should have been calculated at the start of the run, ToDo,
        F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    end
    if ~isnan(CtrlVar.LevelSetDownstreamAGlen)
        F.AGlen(F.LSFMask.NodesOut)=CtrlVar.LevelSetDownstreamAGlen;
        F.n(F.LSFMask.NodesOut)=CtrlVar.LevelSetDownstream_nGlen;
    end
end

if CtrlVar.TestForRealValues
    if ~isreal(F.ub) ; save TestSave ; error('uv:ubNotReal','ub not real!') ; end
    if ~isreal(F.vb) ; save TestSave ; error('uv:vbNotReal','vb not real!') ; end
    if ~isreal(F.ud) ; save TestSave ; error('uv:udNotReal','db not real!') ; end
    if ~isreal(F.vb) ; save TestSave ; error('uv:vdNotReal','vb not real!') ; end
    if ~isreal(l.ubvb) ; save TestSave ; error('uv:ubvbLambdaNotReal','ubvbLambda not real!') ; end
    if ~isreal(l.udvd) ; save TestSave ; error('uv:udvdLambdaNotReal','udvdLambda not real!') ; end
end


if any(F.h<CtrlVar.ThickMin)

    indh0=find(F.h<CtrlVar.ThickMin);
    fprintf('uv: Found ice thicknesses below the min desired ice thickness of hmin=%f, in a diagnostic forward run.\n',CtrlVar.ThickMin)
    fprintf('In total %-i ice nodal values found with too small ice thicknesses, and the min ice thickness in the domain is %f. \n ',numel(indh0),min(F.h));

    if CtrlVar.ResetThicknessToMinThickness==0
        CtrlVar.ResetThicknessToMinThickness=1;
    end


    fprintf('For the purpose of the uv solve, these thickness values will be set to %f \n',CtrlVar.ThickMin)
    [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);

end


Lubvb=[];




switch lower(CtrlVar.FlowApproximation)


    case {'sstream','sstream-rho'}

        if CtrlVar.InfoLevel >= 10 ; fprintf(CtrlVar.fidlog,' Starting SSTREAM diagnostic step. \n') ;  end

        [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]=SSTREAM2dNR2(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
       

        if ~RunInfo.Forward.uvConverged
            fprintf('uv forward calculation did not converge. Resetting ub and vb and solving again.\n')
            F.ub=F.ub*0 ; F.vb=F.vb*0 ; l.ubvb=l.ubvb*0 ;
            % [UserVar,F,l,Kuv,Ruv,RunInfo,Lubvb]=SSTREAM2dNR2(UserVar,CtrlVar,MUA,BCs,F,l,RunInfo);
            [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]=SSTREAM2dNR2(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        end

    case 'ssheet'

        if CtrlVar.InfoLevel >= 10 ; fprintf(CtrlVar.fidlog,' start SSHEET diagnostic. \n') ;  end
        [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);

        [F.ud,F.vd,F.ub,F.vb]=uvSSHEET(CtrlVar,MUA,BCs,F.AGlen,F.n,F.C,F.m,F.rho,F.g,F.s,F.h);
        l.ubvb=[] ; Kuv=[] ; Ruv=[];
        RunInfo.Forward.uvConverged=1;
        RunInfo.Forward.Iterations=NaN;
        RunInfo.Forward.Residual=NaN;

    case 'hybrid'

        if CtrlVar.InfoLevel >= 10 ; fprintf(CtrlVar.fidlog,'Start hybrid: 1:SSTREAM-Step \n') ;  end
        
        [UserVar,F,l,Kuv,Ruv,RunInfo,Lubvb]=SSTREAM2dNR(UserVar,CtrlVar,MUA,BCs,F,l,RunInfo);


        if CtrlVar.InfoLevel >= 10 ; fprintf(CtrlVar.fidlog,' 2:Basal stress, ') ;  end
        GF=GL2d(F.B,F.S,F.h,F.rhow,F.rho,MUA.connectivity,CtrlVar);
        [txzb,tyzb]=CalcNodalStrainRatesAndStresses(CtrlVar,MUA,F.AGlen,F.n,F.C,F.m,GF,F.s,F.b,F.ub,F.vb);
        if CtrlVar.InfoLevel >= 10 ; fprintf(CtrlVar.fidlog,' 3:SSHEET.') ;  end
        [F.ud,F.vd]=uvSSHEETplus(CtrlVar,MUA,BCs,F.AGlen,F.n,F.h,txzb,tyzb);
        if CtrlVar.InfoLevel >= 10 ; fprintf(CtrlVar.fidlog,' Hybrid done \n') ;  end


    case "uvhprescribed"

        [UserVar,RunInfo,F,l]=uvhPrescibed(UserVar,RunInfo,CtrlVar,MUA,F,F,l,BCs);
        Kuv=[] ; Ruv=[] ;
    otherwise

        fprintf('CtrlVar.FlowApproximation %s \n',CtrlVar.FlowApproximation)
        error('what case')


end


% F.h=hOnInput; % this is questionable because then the solve will be for a different h thick

if CtrlVar.InfoLevelCPU> 10
    tdiagnostic=toc(tdiagnostic);
    fprintf(CtrlVar.fidlog,' Ended diagnostic in %-f sec \n ',tdiagnostic) ;
end


if CtrlVar.TestAdjointFiniteDifferenceType=="complex step differentiation"
    CtrlVar.TestForRealValues=false;
end


if  CtrlVar.TestForRealValues
    if ~isreal(F.ub) ; save TestSave ; error('uv:ubNotReal','ub not real!') ; end
    if ~isreal(F.vb) ; save TestSave ; error('uv:vbNotReal','vb not real!') ; end
    if ~isreal(F.ud) ; save TestSave ; error('uv:udNotReal','ud not real!') ; end
    if ~isreal(F.vd) ; save TestSave ; error('uv:vdNotReal','vd not real!') ; end
    if ~isreal(l.ubvb) ; save TestSave ; error('uv:ubvbLambdaNotReal','ubvbLambda not real!') ; end
    if ~isreal(l.udvd) ; save TestSave ; error('uv:udvdLambdaNotReal','udvdLambda not real!') ; end
    if ~isreal(Kuv) ; save TestSave ; error('uv:kvNotReal','kv not real!') ; end
    if ~isreal(Ruv) ; save TestSave ; error('uv:rhNotReal','rh not real!') ; end
end



F.solution="-uv-" ;



if ~RunInfo.Forward.uvConverged

    fprintf("uv solve did not converge. \n  Saving variables in uvDump.mat \n") ;
    save("uvDump.mat","UserVar","RunInfo","CtrlVar","MUA","BCs","F","l")
    F.solution="-no solution-" ;

end






end