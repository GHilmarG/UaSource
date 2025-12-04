function BCsAdjoint=CreatePlausibleBCsForAdjointProblem(BCs,BCsAdjoint)

% Creates `plausible' boundary conditions for the adjoint problem
%
% If the user has already modified BCsAdjoint, do not change
%
% Otherwise:  Set the BCs of the adjoint problem equal to hat of the forward
% problem and then homogenize Dirichlet boundary conditions.
%
% This ensures that all periodic boundary conditions of the forward problem are
% carried over to the adjoint problem. Also the essential boundary conditions of
% both problems are applied over same sections of the computational boundary. 
%
%

BCsDefault=BoundaryConditions;

if isequal(BCsAdjoint,BCsDefault)  % user had not modified BSsAdjoint
    BCsAdjoint=BCs;                % copy BCs of forward run, and the homogenize
    BCsAdjoint.ubFixedValue=BCsAdjoint.ubFixedValue*0;
    BCsAdjoint.vbFixedValue=BCsAdjoint.vbFixedValue*0;
    
end


end