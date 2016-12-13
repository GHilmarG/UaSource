function dJ = CalcBruteForceGradient(func,p0,iRange)



if nargin<3 || isempty(iRange)
    iRange=1:numel(p0);
end


fprintf(' Calculating gradients using brute force method \n')


J0=func(p0);

% Testing gradient using brute force method

delta=1e-10*norm(p0);
dJ=p0*0;

parfor I=iRange
    
    p1=p0;
    p1(I)=p1(I)+delta;
    J1=func(p1);
    dJ(I)=(J1-J0)/delta;
    
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


