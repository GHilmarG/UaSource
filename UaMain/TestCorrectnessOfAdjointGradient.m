


function [J,dJdp,dJdpTest] = TestCorrectnessOfAdjointGradient(func,p0,MUA,CtrlVar,plb,pub,G)



% Get the gradient using the adjoint method
%%
[isA,isB,isC] = isABC(CtrlVar);

[J,dJdp]=func(p0);

NA=MUA.Nnodes;

% Find the subset (iRange) in p, for which the brute-force gradient is to be calculated
if isempty(CtrlVar.Inverse.TestAdjoint.iRange)
    
    stride=10 ; % every 10th node, i.e. 10%
    iRange=(1:stride:MUA.Nnodes)';
else
    iRange=CtrlVar.Inverse.TestAdjoint.iRange;
end

I=(iRange>=1) & (iRange <= MUA.Nnodes);  % Just in case the use sets some CtrlVar.Inverse.TestAdjoint.iRange outside the nodal values in Mesh
iRange=iRange(I);

% if the inversion is done for more than one field, then expand iRange accordingly.
nBlocks=isA+isB+isC;

switch nBlocks
    case 2
        iRange=[iRange(:);iRange(:)+NA];
    case 3
        iRange=[iRange(:);iRange(:)+NA;iRange(:)+2*NA];
end

% select a reasonable step size. This could definitely be improved, but the key thing is to use a constant for log and
% different step size of B compared to a and C.  Here is is also good to try different step sizes for convergence. 
deltaA=0.0001;
deltaB=0.01;
deltaC=0.001; 
if isA
    deltaStepA=deltaA+zeros(MUA.Nnodes,1);
else
    deltaStepA=[];
end
if isB
    deltaStepB=deltaB+zeros(MUA.Nnodes,1) ;
else
    deltaStepB=[];
end
if isC
    deltaStepC=deltaC+zeros(MUA.Nnodes,1);
else
    deltaStepC=[];
end

deltaStep=[deltaStepA;deltaStepB;deltaStepC];


% Gradient calculated using a brute-force finite difference approach
dJdpTest = CalcBruteForceGradient(func,p0,plb,pub,CtrlVar,iRange,deltaStep);

Diff=norm(dJdp(iRange)-dJdpTest(iRange))/norm(dJdp(iRange)+eps);
fprintf("Test Adjoint gradients: Normalized differences between adjoint gradient and FD: %g \n ",Diff)

fig_dJdpTest_bubble=FindOrCreateFigure("Test dJdp") ; clf(fig_dJdpTest_bubble)
plot(dJdp(iRange),dJdpTest(iRange),"or") ; axis equal ;
hold on ;
plot([min(dJdp(iRange)) max(dJdp(iRange))],[min(dJdp(iRange)) max(dJdp(iRange))],"--k")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off
xlabel("Adjoint $\partial J/\partial p$ ",Interpreter="latex")  ;
ylabel("Finite difference $\partial J/\partial p$",Interpreter="latex")
title("$\partial J/\partial p$ :  "+CtrlVar.Inverse.InvertFor,Interpreter="latex")
subtitle(sprintf("Normalized diff %g",Diff),Interpreter="latex")

fig_dJdpTest_bubble=FindOrCreateFigure("Test dJdp bubble") ; clf(fig_dJdpTest_bubble)

bc=bubblechart(dJdp(iRange),dJdpTest(iRange),G(iRange)+1,"r",MarkerFaceColor="g") ; 
axis equal ; 
bubblelegend("$\mathcal{G}+1$",Interpreter="latex",Location="northwest") ;
hold on ; 
plot([min(dJdp(iRange)) max(dJdp(iRange))],[min(dJdp(iRange)) max(dJdp(iRange))],"--k")
xlabel("Adjoint $\partial J/\partial p$ ",Interpreter="latex")  ;
ylabel("Finite difference $\partial J/\partial p$",Interpreter="latex")

drawnow

fprintf("  node                   G        dJdp_adjoint                  dJdp_FD            norm diff        \n")
for iTestNode=1:numel(iRange)
    node=iRange(iTestNode);
    diff=(dJdp(node)-dJdpTest(node))/(abs(dJdp(node)+eps)); 
    fprintf("%7i \t %10.5g \t %15.10g \t %15.10g \t %15.10g \n",node,G(node),dJdp(node),dJdpTest(node),diff)
end

%%

end
%%