   
function PlotCauchy2NewtonPath(CtrlVar,xVector,fVector,R2C,R2N,gammaminCN,rCNmin,R20)

xVector=[xVector;0;1]; fVector=[fVector;R2C;R2N];
I=~isnan(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;
[~,I]=sort(xVector)  ; xVector=xVector(I) ; fVector=fVector(I) ;


fig=FindOrCreateFigure("Cauchy2NewtonPath") ; clf(fig)
plot(xVector,fVector,'o-r') ;
hold on ;
plot(0,R2C,"*m")
plot(1,R2N,"*b")


text(0,R2C,"   C",HorizontalAlignment="left")
text(1,R2N,"N   ",HorizontalAlignment="right")

yline(R20,"--k","$r_0$",interpreter="latex",LabelHorizontalAlignment="center")

%yline(rminNewton,"-.k","$\min r_N$",interpreter="latex",LabelHorizontalAlignment="right")
%yline(r2MD,"-.k","$\min r_C$",interpreter="latex",LabelHorizontalAlignment="left")



plot(gammaminCN,rCNmin ...
    ,'o',MarkerFaceColor="b",MarkerSize=10)

xlabel("Distance along the Cauchy-to-Newton path",interpreter="latex")
ylabel("$r^2$",interpreter="latex")

legend("$\\|R\\|^2$","$\\|R_C\\|^2$","$\\|R_N\\|$","$r^2_0$","$\min r^2_{CN}$",interpreter="latex",location="best")
title("From the Cauchy-Point to the Newton-Point")

if isfield(CtrlVar,"BacktrackIteration")
    subtitle(sprintf("Iteration:%i",CtrlVar.BacktrackIteration))
end



if isfield(CtrlVar,"time")
    subtitle(sprintf("t=%f   dt=%f",CtrlVar.time,CtrlVar.dt),Interpreter="latex")
end
% fig=gcf ; exportgraphics(fig,"ExampleCaucyToNewtonPath.pdf")
drawnow


end