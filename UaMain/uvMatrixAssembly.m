



function [RunInfo,Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F,BCs)

%
% Ruv=Tint-Fext;
% Tint   : internal nodal forces
% Fint   : external nodal forces

narginchk(4,5)
nargoutchk(2,4)



if nargout==2 && ~CtrlVar.uvMatrixAssembly.Ronly
    error("uvMatrixAssembly: More than one output, but assembly only required for R. ")
end

if nargin==4
    BCs=[];
end


poolobj = gcp('nocreate');  % check if parpool exists, but do not create one if it does not exist already

if isempty(poolobj)

    if CtrlVar.Parallel.uvAssembly.spmd.isOn

        
       % fprintf("\n ======= No parallel pool is open. To run %ca using parallel options, a parallel pool must be opened ahead of a call to %ca.\n",218,218)
       % fprintf(" ======= Parallel options are turned off.\n")

        CtrlVar.Parallel.uvhAssembly.parfor.isOn=false;
        CtrlVar.Parallel.uvhAssembly.spmd.isOn=false;
        CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=[];

        CtrlVar.Parallel.uvAssembly.spmd.isOn=false;
        CtrlVar.Parallel.uvAssembly.parfeval.isOn=false;

        CtrlVar.Parallel.uvAssembly.spmd.nWorkers=[];

        CtrlVar.Parallel.isTest=false;

        CtrlVar.Parallel.hAssembly.parfor.isOn=false ;
        CtrlVar.Parallel.LSFAssembly.parfor.isOn=0;

        CtrlVar.Parallel.Distribute=false;

    end

end




tAssembly=tic;

switch lower(CtrlVar.FlowApproximation)

    case "sstream"


        if CtrlVar.Parallel.uvAssembly.spmd.isOn


            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,BCs);

        else  % this is the otherwise default sequential assembly

            [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs) ;

        end



    case "sstream-rho"

        [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAMrho(CtrlVar,MUA,F,BCs) ;

    otherwise

        error("UA:uvMatrixAssemblyCaseNotFound","Case Not Found")

end


RunInfo.CPU.Assembly.uv=RunInfo.CPU.Assembly.uv+toc(tAssembly);

end



