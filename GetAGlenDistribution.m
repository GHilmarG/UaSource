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



if CtrlVar.LevelSetMethod &&  ~isnan(CtrlVar.LevelSetDownstreamAGlen )

    if isempty(F.LSFMask)
        F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    end
    F.AGlen(F.LSFMask.NodesOut)=CtrlVar.LevelSetDownstreamAGlen;


    % F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    % F.AGlen(F.LSFMask.NodesOut)=CtrlVar.LevelSetDownstreamAGlen;

end





end