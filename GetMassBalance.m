function [UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)


nOut=nargout;
if nOut~=5
    error('Ua:GetMassBalance','Need 5 output arguments')
end


if CtrlVar.MassBalanceGeometryFeedback>0
    
    N=nargout('DefineMassBalance');
    
    switch N
        
        case 4
            
            [as,ab,dasdh,dabdh]=DefineMassBalance(CtrlVar.Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
            
        case 5
            
            [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
            
    end
    
else
    
    N=nargout('DefineMassBalance');
    
    switch N
        
        case 2
            [as,ab]=DefineMassBalance(CtrlVar.Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
            
        case 3
            
            [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
            
    end
    
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