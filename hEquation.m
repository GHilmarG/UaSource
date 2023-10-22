





function [UserVar,h,lambda]=hEquation(UserVar,CtrlVar,MUA,F,BCs,Priors,Meas)

%%
% Solves the linear h-equation on the form:
%
% $$ d (u h)/dx + d (v h)/dy - \nabla \cdot (\kappa \nabla h) = a$$
%
% for h
%
% The natural boundary condition is
%
%  $$\nabla h \cdot \hat{n} = 0 $$
%
% ie,  the free outflow condition
%
% Boundary conditions: The BCs are identical to define thickness boundary
% conditions.
%
%
%
%
% * Method="Fh +  <l , hMeas>";
%
% Forward model solved for measurements as boundary conditions. This can be thought of as using both the model and the
% thickness measurements as hard constraints.
%
%
%
% * Method="(h-hmeas) P (h-hmeas) / 2 +  <l , Fh>"
%
% Deviation from measurments minimized with forward model as a constraint, i.e. measurments are a soft constraint and the
% model a hard constraint.
%
% system to solve
%
%  [P K^T ] [h] = [P hmeas]
%  [K   0 ] [l]   [ b]
%
% where K= D_h Fh and the linearized forward problem therefore on the form K h = b
%
%%

narginchk(7,7)


CtrlVar.hEq.Method=replace(CtrlVar.hEq.Method,"  "," ");  % get rid of any double spaces in the string Method 

switch CtrlVar.hEq.Method

    case "Fh + <l , hMeas>"

        % Forward model, Fh=0, and measurements are hard constraints

        [UserVar,K,b]=hEquationAssembly(UserVar,CtrlVar,MUA,F.ub,F.vb,F.as);


        MLC=BCs2MLC(CtrlVar,MUA,BCs);
        L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;


        h0=F.x*0;
        [h,lambda]=solveKApe(K,L,b,Lrhs,h0,lambda,CtrlVar);
        h=full(h);

    case "(h-hmeas) P (h-hmeas) / 2 + <l , Fh>"


        % Measurements soft constraints
        % Model a soft constraint

        % for the time being the errors are prescribed here directly
        Error=F.s*0+1e10;        % Errors where no measurements are avalable set to a high value of 1e10
        hmeas=F.s*0;             % Simply set thickness 'measurments' to the value of 0 where none are available

        Error(BCs.hFixedNode)=1;                 % Errors where measurements are avalable set to a value of 1
        hmeas(BCs.hFixedNode)=BCs.hFixedValue;   % Here the BCs are introduced as a soft contraint


        P=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1./(Error.^2),MUA.Nnodes,MUA.Nnodes);


        % system to solve
        %
        %  [P K^T ] [h] = [P hmeas]
        %  [K   0 ] [l]   [ b]
        %

        c=P*hmeas ;

        [UserVar,K,b]=hEquationAssembly(UserVar,CtrlVar,MUA,F.ub,F.vb,F.as);


        h0=F.x*0; lambda=F.x*0 ;
        [h,lambda]=solveKApe(P,K,c,b,h0,lambda,CtrlVar);
        h=full(h);

    case "(h-hmeas) P (h-hmeas) / 2 + (h-hprior) Q (h-hprior) / 2 + (K h - b) M (K h - b) / 2"

        % h measurements precision matrix: P 
        hMeasPrecision=1./Meas.hCov; 
        [UserVar,P]=GeneralizedMassAndStiffnessAssembly(UserVar,CtrlVar,MUA,hMeasPrecision,0,0,0);


        % h prior precision matrix: Q
        if ~isfield(MUA,"Dxx")
            [MUA.Dxx,MUA.Dyy]=StiffnessMatrix2D1dof(MUA);
        end
        

        gha=CtrlVar.hEq.gha ;
        ghs=CtrlVar.hEq.ghs ;
        Q=gha* MUA.M+ghs*(MUA.Dxx+MUA.Dyy) ;



        % Forward model  precision matrix: M
        gFa=CtrlVar.hEq.gFa ; 
        M=gFa*MUA.M ;

        % forward model assembly
        [UserVar,K,b]=hEquationAssembly(UserVar,CtrlVar,MUA,F.ub,F.vb,F.as);

        A=P+Q+K'*M*K ;
        B=P*Meas.h+Q*Priors.h+K'*M*b ;

        h=A\B ;
        h=full(h);
        lambda=[] ; 


    otherwise

        error("case not found")


end




