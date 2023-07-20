function [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F)

narginchk(4,4)
nargoutchk(2,2)


InputFile="DefineMassBalance.m" ; TestIfInputFileInWorkingDirectory(InputFile) ;



N=nargout('DefineMassBalance');
NargInputFile=nargin(InputFile);

switch N
    
    case 3
        
        
        if NargInputFile>4
            
            [UserVar,F.as,F.ab]=DefineMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
            
        else
            
            [UserVar,F.as,F.ab]=DefineMassBalance(UserVar,CtrlVar,MUA,F);
            
        end
        
        
        F.dasdh=F.as*0 ;  F.dabdh=F.ab*0 ;
        
    case 5
        
        if NargInputFile>4
            
            [UserVar,F.as,F.ab,F.dasdh,F.dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
            
        else
            
            [UserVar,F.as,F.ab,F.dasdh,F.dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,F);
            
        end
        
        
    otherwise
        
        fprintf('DefineMassBalance must return either 3 or 5 outputs\n')
        fprintf('The outputs must be either: \n')
        fprintf('\t (UserVar,as,ab) \n')
        fprintf('or  \n')
        fprintf('\t (UserVar,as,ab,dasdh,dabdh) \n')
        error('Ua:IncorrectUserInputs','Incorrect number of outputs returned by DefineMassbalance.m')
        
end


% some input checks

errorStruct.identifier = 'GetMassBalance:NaNinInput';
if any(isnan(F.as))
    errorStruct.message = 'nan in as';
    error(errorStruct)
end

if any(isnan(F.ab))
    errorStruct.message = 'nan in ab';
    error(errorStruct)
end

if any(isnan(F.dasdh))
    errorStruct.message = 'nan in dasdh';
    error(errorStruct)
end

if any(isnan(F.dabdh))
    errorStruct.message = 'nan in dabdh';
    error(errorStruct)
end

if numel(F.as)==1
    F.as=F.as+zeros(MUA.Nnodes,1);
end

if numel(F.ab)==1
    F.ab=F.ab+zeros(MUA.Nnodes,1);
end

if numel(F.dasdh)==1
    F.dasdh=F.dasdh+zeros(MUA.Nnodes,1);
end

if numel(F.dabdh)==1
    F.dabdh=F.dabdh+zeros(MUA.Nnodes,1);
end


if  MUA.Nnodes ~= numel(F.as)
    fprintf('as must have same number of values as there are nodes in the mesh \n')
    error('DefineMassBalance returns incorrect dimensions ')
end


if  MUA.Nnodes ~= numel(F.ab)
    fprintf('ab must have same number of values as there are nodes in the mesh \n')
    error('DefineMassBalance returns incorrect dimensions ')
end

if ~iscolumn(F.ab)
    error("GetMassBalance:abNotColumnVector","error: ab returned by DefineMassBalance not a column vector.")
end

if ~iscolumn(F.as)
    error("GetMassBalance:abNotColumnVector","error: as returned by DefineMassBalance not a column vector.")
end




end