






function  [UserVar,ub,vb,ud,vd]=DefineStartVelValues(UserVar,CtrlVar,MUA,BCs,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m)
          

%%
% Define start values for velocities
%
% [UserVar,ub,vb,ud,vd]=DefineStartVelValues(UserVar,CtrlVar,MUA,BCs,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m)
%
% This user m-file defines the starting values for the velocities.  This is just
% the initial estimate of the solution and, provided the solver converges, it has
% no impact on the final solution. On the other hand, if a good initial estimate is
% available, then prescribing it may speed up the solution. In most cases a good
% initial estimate is simply to set all velocities to zero and this is the
% default approach. So generally this m-file is not required to obtain a solution,
% but it may speed things up.
%
%


        
end
