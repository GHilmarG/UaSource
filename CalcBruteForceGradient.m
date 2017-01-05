function dJ = CalcBruteForceGradient(func,p0,CtrlVar)



if isempty(CtrlVar.Inverse.TestAdjoint.iRange)
    iRange=1:numel(p0);
else
    iRange=CtrlVar.Inverse.TestAdjoint.iRange;
end


fprintf(' Calculating gradients using brute force method. \n')


J0=func(p0);

% Testing gradient using brute force method

delta=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize*norm(p0);

dJ=p0*0+NaN;


switch  lower(CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType)
    
    case 'forward'
        
        parfor I=iRange
            
            p1=p0;
            p1(I)=p1(I)+delta;
            J1=func(p1);
            dJ(I)=(J1-J0)/delta;
            
        end
        
    case 'central'
        
        parfor I=iRange
            p1=p0;
            pm1=p0;
            p1(I)=p1(I)+delta;
            J1=func(p1);
            pm1(I)=pm1(I)-delta;
            Jm1=func(pm1);
            
            dJ(I)=(J1-Jm1)/delta/2;
        end
    otherwise
        
        fprintf(' CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType has invalid value. \n')
        error('which case')
end

% %%
% calculates d kv/ d b using the complex number method
% using the fact that df/dx=Im(f(x+i dx))/dx
% delta=1e-6*norm(p0);
% dJ=p0*0;
% 
% for I=iRange
%     
%     p1=p0;
%     p1(I)=p1(I)+1i*delta;
%     J1=func(p1);
%     dJ(I)=imag(J1)/delta;
%     
% end

%%
%DkvDb=imag(kv)/db;
%DrhDb=imag(rh)/db;

end


