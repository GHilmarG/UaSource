
function [slope0,gammaMin,Q]=QgHlsq(F,J,Aeq,beq,x0,l0,dx,dl,gamma)

%%
%
% The slopes at zero and gamma min are same for both the constrained and the unconstrained quadratic problems
%
%
%%



arguments
    F (:,1) {mustBeNumeric}
    J (:,:) {mustBeNumeric}
    Aeq (:,:) {mustBeNumeric}
    beq (:,1) {mustBeNumeric}
    x0 (:,1) {mustBeNumeric}
    l0 (:,1) {mustBeNumeric}
    dx (:,1) {mustBeNumeric}
    dl (:,1) {mustBeNumeric}
    gamma (1,1) {mustBeNumeric}
end

if isempty(Aeq)

    % Note that the search direction is [dx;dl]
    slope0=2*F'*J*dx;

else
    slope0=2*F'*J*dx+(Aeq'*l0)'*dx+(Aeq*x0 -beq)'*dl;
end

slope0=full(slope0);
gammaMin=-slope0/((J*dx)'*(J*dx));


if ~isnan(gamma)


    dx=gamma*dx;
    dl=gamma*dl;

    Q=0.5*(J*dx)'*(J*dx)+F'*(J*dx) ;

    if ~isempty(Aeq)
        Q=Q +0.5*(dx'*Aeq'*dl+dl'*Aeq*dx) + (Aeq'*l0)'*dx+(Aeq*x0-beq)'*dl;
    end

else
    Q=nan;
end





end
