function [as,ab,dasdh,dabdh]=GetMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)


if CtrlVar.MassBalanceGeometryFeedback>0
    [as,ab,dasdh,dabdh]=DefineMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
else
    [as,ab]=DefineMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
    dasdh=as*0 ;  dabdh=ab*0 ;
end


% some input checks

errorStruct.identifier = 'GetMassBalance:NaNinInput';
if any(isnan(as))
    errorStruct.message = 'nan in as';
    error(errorStruct)
end

if any(isnan(ab))
    errorStruct.message = 'nan in ab';
    error(errorStruct)
end

if any(isnan(dasdh))
    errorStruct.message = 'nan in dasdh';
    error(errorStruct)
end

if any(isnan(dabdh))
    errorStruct.message = 'nan in dabdh';
    error(errorStruct)
end


if numel(as)==1
    as=as+zeros(MUA.Nnodes,1);
end

if numel(ab)==1
    ab=ab+zeros(MUA.Nnodes,1);
end


if  MUA.Nnodes ~= numel(as)
    fprintf('as must have same number of values as there are nodes in the mesh \n')
    error('DefineMassBalance returns incorrect dimentions ')
end


if  MUA.Nnodes ~= numel(ab)
    fprintf('ab must have same number of values as there are nodes in the mesh \n')
    error('DefineMassBalance returns incorrect dimentions ')
end



end