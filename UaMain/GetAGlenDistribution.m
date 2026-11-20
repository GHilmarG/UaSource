function [UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F)


narginchk(4,4)
nargoutchk(2,2)

InputFile="DefineAGlenDistribution.m"; TestIfInputFileInWorkingDirectory(InputFile) ;
NargInputFile=nargin(InputFile);

N=nargout('DefineAGlenDistribution');


switch N
    
    case 2
        
        
        if NargInputFile>4
            [F.AGlen,F.n]=DefineAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
        else
            [F.AGlen,F.n]=DefineAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUA,F);
        end
        
    case 3
        
        if NargInputFile>4
            [UserVar,F.AGlen,F.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
        else
            [UserVar,F.AGlen,F.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,F);
        end
        
    otherwise
        
        error('Ua:GetAGlen','DefineAGlenDistribution must return either 2 or 3 output arguments')
        
end

[F.AGlen,F.n]=TestAGlenInputValues(CtrlVar,MUA,F.AGlen,F.n);


F.AGlenmax=CtrlVar.AGlenmax;
F.AGlenmin=CtrlVar.AGlenmin;


%% Better do these test in uvh and uv, when they are needed


% (20 Nov 2026) Values outside [AGlenmin,AGlenmax] are reported at CtrlVar.InfoLevel>=1 (ReportClipping.m). Previously the report
% was only printed at InfoLevel>=10, and with incorrect formatting for values above AGlenmax.
ReportClipping(CtrlVar,"AGlen",F.AGlen,F.AGlenmin,F.AGlenmax,LowerName="CtrlVar.AGlenmin",UpperName="CtrlVar.AGlenmax",Where="GetAGlenDistribution") ;
F.AGlen=kk_proj(F.AGlen,F.AGlenmax,F.AGlenmin);



if CtrlVar.LevelSetMethod &&  ~isnan(CtrlVar.LevelSetDownstreamAGlen) &&  ~isnan(CtrlVar.LevelSetDownstream_nGlen)
    if isempty(F.LSFMask)  % This should have been calculated at the start of the run, ToDo,
        F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    end
    if ~isnan(CtrlVar.LevelSetDownstreamAGlen)
        F.AGlen(F.LSFMask.NodesOut)=CtrlVar.LevelSetDownstreamAGlen;
        F.n(F.LSFMask.NodesOut)=CtrlVar.LevelSetDownstream_nGlen;
    end
end

  




end