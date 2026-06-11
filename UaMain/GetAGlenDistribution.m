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


[F.AGlen,iU,iL]=kk_proj(F.AGlen,F.AGlenmax,F.AGlenmin);

niU=numel(find(iU));
niL=numel(find(iL));

if CtrlVar.InfoLevel>=10
    if niU>0
        fprintf('TestAGlenInputValues: On input %i AGlen values are larger than largest allowed value (%g).  \n',niU);
        fprintf(' These values have been reset to the maximum allowed value of %g .\n',niU,CtrlVar.AGlenmax);

    end
    if niL>0
        fprintf('TestAGlenInputValues: On input %i AGlen values are smaller than smallest allowed value (%g). \n',niL,CtrlVar.AGlenmin);
        fprintf(' These values are have been reset to the minimum allowed value of %g.\n',CtrlVar.AGlenmin);
    end
end



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