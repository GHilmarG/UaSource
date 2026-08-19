
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
    slope0=F'*J*dx;
    Jdx=J*dx; 
    gammaMin=-slope0/((Jdx)'*(Jdx));

else

    % If the iterate is feasible the slopes at origin, and the step size to the minimum, are identical for the constraint and
    % unconstrained problems. However for numerical reasons I use the full expression here.

    slope0=F'*J*dx+(Aeq'*l0)'*dx+(Aeq*x0 -beq)'*dl;  
    Jdx=J*dx; 
    gammaMin=-slope0/((Jdx)'*(Jdx)+2*(Aeq*dx)'*dl);

end

slope0=full(slope0);

Q=Qlsq(F,J,Aeq,beq,x0,l0,dx,dl,gamma);



% doPlots=true;
% 
% % Plot Q along [dx;dl] direction
% 
% gammaMax=1;
% gammaVector=linspace(0,gammaMax,100); gammaVector=gammaVector(:); 
% Qvector=nan(numel(gammaVector),1);
% 
% 
% for I=1:numel(gammaVector)
%     Qvector(I)=Qlsq(F,J,Aeq,beq,x0,l0,dx,dl,gammaVector(I)) ;
% end
% FigQ=FindOrCreateFigure("Qlsq along search direcion") ; clf(FigQ)
% 
% plot(gammaVector,Qvector)
% 
% 
% 


end





function Q=Qlsq(F,J,Aeq,beq,x0,l0,dx,dl,g)

if ~isnan(g)


    dx=g*dx;
    dl=g*dl;

    Q=0.5*(J*dx)'*(J*dx)+F'*(J*dx) ;

    if ~isempty(Aeq)
        Q=Q +0.5*(dx'*Aeq'*dl+dl'*Aeq*dx) + (Aeq'*l0)'*dx+(Aeq*x0-beq)'*dl;
    end

else
    Q=nan;
end


end





