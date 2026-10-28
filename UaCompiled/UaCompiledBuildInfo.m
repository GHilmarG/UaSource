function [ok,Reason,Note]=UaCompiledBuildInfo(Action,MexNames)

%% Build information of the MEX files in UaCompiled: records it when building, and checks it before use
%
%   UaCompiledBuildInfo("record",MexNames)                   % called by BuildUaCompiled.m after building the MEX files MexNames
%   [ok,Reason,Note]=UaCompiledBuildInfo("check",MexNames)   % called by uvhCompiledKernelCheck.m (once per session and variant)
%
% For each MEX file, the MATLAB release it was built with, and a checksum of each source file compiled into it, are kept in
% UaCompiledBuildInfo_<mexext>.mat in this folder (one file per platform, eg UaCompiledBuildInfo_mexw64.mat).
%
% "check" returns ok=false, with the reason, if
%   - the MEX file was built with another MATLAB release than the one running, or
%   - one of its source files has changed since it was built (the files are named in the reason).
% The vectorised assembly is then used instead (uvhMatrixAssembly.m), until the MEX files are rebuilt with BuildUaCompiled.m.
%
% If there is no build information for a MEX file (eg built before 26 Oct 2026), ok=true and Note says that it could not be checked.
%
% The checksums (ArrayChecksum.m) are of the file contents with carriage returns removed, so that a checkout with other line endings
% (git core.autocrlf) does not count as a change. Source files are recorded relative to the UaSource folder.
%
% See also: BuildUaCompiled, uvhCompiledKernelCheck
%
% (26 Oct 2026)
%%

arguments
    Action (1,1) string {mustBeMember(Action,["record","check"])}
    MexNames string
end

ok=true ; Reason="" ; Note="" ;

Folder=string(fileparts(mfilename('fullpath'))) ;   % UaSource/UaCompiled
Root=string(fileparts(Folder)) ;                    % UaSource
InfoFile=fullfile(Folder,"UaCompiledBuildInfo_"+mexext+".mat") ;

switch Action

    case "record"

        if isfile(InfoFile)
            S=load(InfoFile,"Info") ; Info=S.Info ;
        else
            Info=struct("Mex",{},"Release",{},"Built",{},"Files",{},"Checksums",{}) ;
        end
        for m=MexNames(:)'
            Files=SourceFiles(m,Root) ;
            Checksums=zeros(numel(Files),3) ;
            for k=1:numel(Files)
                Checksums(k,:)=FileChecksum(fullfile(Root,Files(k))) ;
            end
            Entry=struct("Mex",m,"Release",string(version('-release')),"Built",string(datetime('now')),"Files",Files,"Checksums",Checksums) ;
            i=find(string([Info.Mex])==m,1) ;
            if isempty(i) ; Info(end+1)=Entry ; else ; Info(i)=Entry ; end %#ok<AGROW>
        end
        save(InfoFile,"Info") ;

    case "check"

        if ~isfile(InfoFile)
            Note="there is no build information for the MEX files ("+InfoFile+"), so it could not be checked whether they are up to date. Rebuilding them with BuildUaCompiled.m creates it" ;
            return
        end
        S=load(InfoFile,"Info") ; Info=S.Info ;
        Release=string(version('-release')) ;
        Changed=strings(0) ; OtherRelease=strings(0) ; NoInfo=strings(0) ;
        for m=MexNames(:)'
            i=find(string([Info.Mex])==m,1) ;
            if isempty(i) ; NoInfo(end+1)=m ; continue ; end %#ok<AGROW>
            if Info(i).Release~=Release ; OtherRelease(end+1)=m+" (R"+Info(i).Release+")" ; end %#ok<AGROW>
            for k=1:numel(Info(i).Files)
                f=fullfile(Root,Info(i).Files(k)) ;
                if ~isfile(f) || ~isequal(FileChecksum(f),Info(i).Checksums(k,:))
                    Changed(end+1)=Info(i).Files(k) ; %#ok<AGROW>
                end
            end
        end
        Changed=unique(Changed,"stable") ;
        if ~isempty(OtherRelease)
            ok=false ;
            Reason="the MEX files were built with another MATLAB release than R"+Release+": "+strjoin(OtherRelease,", ")+" (rebuild them with BuildUaCompiled.m)" ;
        elseif ~isempty(Changed)
            ok=false ;
            Reason="the MEX files are outdated: "+strjoin(Changed,", ")+" changed since they were built (rebuild them with BuildUaCompiled.m)" ;
        end
        if ~isempty(NoInfo)
            Note="there is no build information for "+strjoin(NoInfo,", ")+", so it could not be checked whether they are up to date. Rebuilding them with BuildUaCompiled.m creates it" ;
        end

end

end

%% local functions

function Files=SourceFiles(MexName,Root)
% Source files compiled into the MEX file MexName, relative to Root (UaSource)
if MexName=="SparseFromCSC"
    Full=string(fullfile(Root,"UaCompiled","SparseFromCSC.c")) ;
else
    % the kernel: the entry points and all UaSource functions they call, and MakeKernelParamsFlex.m (its list of fields fixes the
    % type of the kernel's parameter input). FindOrCreateFigure.m is only called by plotting that is excluded from compiled code.
    Full=string(matlab.codetools.requiredFilesAndProducts({'uvhKernelFull','uvhKernelRonly'}))' ;
    Full=[Full ; string(which('MakeKernelParamsFlex'))] ;
    Full=Full(startsWith(Full,Root+filesep) & ~endsWith(Full,filesep+"FindOrCreateFigure.m")) ;
end
Files=unique(erase(Full,Root+filesep)) ;
Files=replace(Files,"\","/") ;   % the same on all platforms
end

function h=FileChecksum(f)
% Checksum of the contents of file f, with carriage returns removed (line endings do not matter)
fid=fopen(f,'r') ;
b=fread(fid,inf,'*uint8') ;
fclose(fid) ;
b(b==13)=[] ;
h=ArrayChecksum(b) ;
end
