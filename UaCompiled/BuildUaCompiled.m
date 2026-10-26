function BuildUaCompiled(Variants)

%% Builds the MEX files of the compiled uvh assembly (requires MATLAB Coder and a C compiler, preferably with OpenMP)
%
%   BuildUaCompiled                    % builds all MEX files
%   BuildUaCompiled(["Full_n3i6","Sparse"])   % builds only the listed ones
%
% Variants:
%
%   Full, Ronly            : kernel, residual and Jacobian / residual only, any element type with nod<=10 and nip<=37 (bounded sizes)
%   Full_n3i6, Ronly_n3i6  : the same for linear elements with 6 integration points only (fixed sizes, faster)
%   Sparse                 : SparseFromCSC, creates the sparse matrix directly from its CSC arrays (hand-written C, built with mex)
%
% The MEX files are written to this folder, the intermediate files of the code generation to tempdir.
%
% A MEX file that has been called in the current MATLAB session cannot be replaced (the OpenMP threads keep it loaded),
% so build in a fresh session, for example:
%
%   matlab -batch "addpath(genpath('<UaSource>')); BuildUaCompiled"
%
% See README.md in this folder.
%
%%

if nargin<1 || isempty(Variants)
    Variants=["Full","Ronly","Full_n3i6","Ronly_n3i6","Sparse"];
end
Variants=string(Variants);

Folder=fileparts(mfilename('fullpath'));

cfg=coder.config('mex');
cfg.EnableOpenMP=true;
cfg.IntegrityChecks=false;
cfg.ResponsivenessChecks=false;

% Type of the run-time parameter structure (the field set and types are fixed, see MakeKernelParamsFlex.m)
R=MakeKernelParamsFlex(Ua2D_DefaultParameters());
RT=coder.typeof(R);
RT.Fields.SlidingLaw=coder.typeof('a',[1 Inf]);
RT.Fields.ThicknessPenaltyMassBalanceFeedbackFunction=coder.typeof('a',[1 Inf]);
RT.Fields.uvh.Fields.SUPG.Fields.tau=coder.typeof('a',[1 Inf]);

% Mesh structure (without Deriv and DetJ, which are passed as separate arrays), nodal vectors, index arrays
MUAt=coder.typeof(struct('Nele',0,'nod',0,'nip',0,'points',zeros(6,2),'weights',zeros(1,6),'EleAreas',zeros(2,1)));
MUAt.Fields.EleAreas=coder.typeof(0,[Inf 1],[1 0]);
V=coder.typeof(0,[Inf 1],[1 0]);   % nodal vectors, Nnodes x 1 (q, muk, V0 and the melange fields may be 0 x 1)
U=coder.typeof(uint32(0),[Inf 1],[1 0]);   % invPtr, invIdx (see uvhCSCPattern.m)

for v=Variants

    switch v

        case "Sparse"

            fprintf('BuildUaCompiled: SparseFromCSC.c -> SparseFromCSC\n');
            cc=mex.getCompilerConfigurations('C','Selected');
            if isempty(cc) ; error('BuildUaCompiled:NoCompiler','No C compiler selected for mex (see mex -setup C).') ; end
            if contains(cc.ShortName,["mingw","gcc"],'IgnoreCase',true)
                Flags={'CFLAGS=$CFLAGS -fopenmp -O2','LDFLAGS=$LDFLAGS -fopenmp'};
            elseif contains(cc.ShortName,"MSVC",'IgnoreCase',true)
                Flags={'COMPFLAGS=$COMPFLAGS /openmp'};
            else
                Flags={};   % without OpenMP: the arrays are then copied serially
                fprintf('BuildUaCompiled: compiler %s, SparseFromCSC built without OpenMP\n',cc.ShortName);
            end
            mex('-R2018a','-silent',Flags{:},'-outdir',Folder,fullfile(Folder,'SparseFromCSC.c'));

        otherwise   % kernel variants

            if endsWith(v,"_n3i6")
                nod=3; nip=6; vs=[0 0];      % fixed sizes
            else
                nod=10; nip=37; vs=[1 1];    % upper bounds
            end
            M=MUAt;
            M.Fields.points=coder.typeof(0,[nip 2],[vs(2) 0]);
            M.Fields.weights=coder.typeof(0,[1 nip],[0 vs(2)]);
            Deriv=coder.typeof(0,[Inf 2 nod nip],[1 0 vs(1) vs(2)]);
            DetJ=coder.typeof(0,[Inf nip],[1 vs(2)]);
            Conn=coder.typeof(0,[Inf nod],[1 vs(1)]);
            % MUA, Deriv, DetJ, connectivity, b h u v AGlen n C m (8), q muk V0 (3), h0 u0 v0 as0 ab0 as1 ab1 dadh B S rho (11), LSFMask hBCsMask (2),
            % melange (8), R, rhow g ca sa dt, IntOut, invPtr, invIdx
            args=[{M,Deriv,DetJ,Conn}, repmat({V},1,8), {V,V,V}, repmat({V},1,11), {V,V}, repmat({V},1,8), {RT,0,0,0,0,0,false,U,U}];
            Entry="uvhKernel"+extractBefore(v+"_","_");
            OutName="uvhKernel"+v+"_mex";
            fprintf('BuildUaCompiled: %s -> %s\n',Entry,OutName);
            codegen(char(Entry),'-config',cfg,'-args',args,'-o',char(fullfile(Folder,OutName)),'-d',char(fullfile(tempdir,'UaCompiledBuild',OutName)));

    end
end

end
