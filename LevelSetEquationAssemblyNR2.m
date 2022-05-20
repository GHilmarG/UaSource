function [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,F0,F1)

% [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)


narginchk(5,5)
nargoutchk(2,6)

if ~(nargout==2 || nargout==6)
    error('LevelSetEquationAssemblyNR2:IncorrectNumberOfInputs','Incorrect number of input arguments')
end

[UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,F0,F1);

% f1=F1.LSF ;  f0=F0.LSF ; 
% [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);





end