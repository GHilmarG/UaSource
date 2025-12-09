



          
function [UserVar,ub,vb,ud,vd,l]=DefineStartVelValues(UserVar,CtrlVar,MUA,BCs,F,l) 

%%
% Define start values for velocities
%
%  [UserVar,F.ub,F.vb,F.ud,F.vd]=DefineStartVelValues(UserVar,CtrlVar,MUA,BCs,F) 
%
% This user m-file defines the starting values for the velocities.  This is just
% the initial estimate of the solution and, provided the solver converges, it has
% no impact on the final solution. On the other hand, if a good initial estimate is
% available, then prescribing it may speed up the solution. In most cases a good
% initial estimate is simply to set all velocities to zero and this is the
% default approach. So generally this m-file is not required to obtain a solution,
% but it may speed things up.
%
%  l are Lagrange parameters related to boundary conditions
%
% If you have those from a previous solve, these can be specified here.
%
%%

ub=F.ub;
vb=F.vb;

ud=F.ud;
vd=F.vd;

        
end
