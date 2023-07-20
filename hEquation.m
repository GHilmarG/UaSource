function [UserVar,h,lambda]=hEquation(UserVar,CtrlVar,MUA,F,BCs,kIso,kAlong,kCross,Method)

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

if nargin< 9
    Method="Fh +  <l hMeas>";
end


switch Method

    case "Fh +  <l hMeas>"


        [UserVar,K,b]=hEquationAssembly(UserVar,CtrlVar,MUA,F.ub,F.vb,F.as,kIso,kAlong,kCross);

 
        MLC=BCs2MLC(CtrlVar,MUA,BCs);
        L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;


        h0=F.x*0;
        [h,lambda]=solveKApe(K,L,b,Lrhs,h0,lambda,CtrlVar);
        h=full(h);

    case "(h-hmeas) P (h-hmeas) / 2 +  <l , Fh>"


        Error=F.s*0+1e10; 
        hmeas=F.s*0;

        Error(BCs.hFixedNode)=1;
        hmeas(BCs.hFixedNode)=BCs.hFixedValue;
        

        P=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1./(Error.^2),MUA.Nnodes,MUA.Nnodes); 

        % system to solve
        %
        %  [P K^T ] [h] = [P hmeas]
        %  [K   0 ] [l]   [ b]
        %

        c=P*hmeas ;

         [UserVar,K,b]=hEquationAssembly(UserVar,CtrlVar,MUA,F.ub,F.vb,F.as,kIso,kAlong,kCross);


        h0=F.x*0; lambda=F.x*0 ; 
        [h,lambda]=solveKApe(P,K,c,b,h0,lambda,CtrlVar);
        h=full(h);




end




