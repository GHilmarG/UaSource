function PlotForceResidualVectors2(CtrlVar,MUA,F,msg,R,L,lambda,iteration)


x=MUA.coordinates(:,1);
y=MUA.coordinates(:,2);

Nnodes=length(MUA.coordinates);
if ~isempty(L)
    R=R+L'*lambda;
end

FigName='Nodal force residuals' ;
figNRr=FindOrCreateFigure(FigName); clf(figNRr) ;


if ~contains(msg,'h-only')  % uvh residuals
    % uv-residuals

    quiver(x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,R(1:Nnodes),R(Nnodes+1:2*Nnodes))
    title("Nodal-force residuals (R+L^T \lambda)")

    % h residuals
    if contains(msg,'h')
        hold on

        Rmax=max(abs(R(2*Nnodes+1:end))) ;
        L=min([max(x)-min(x),max(y)-min(y)])/CtrlVar.PlotXYscale ;
        alpha=0.01 ;
        PlotScale=Rmax/(alpha*L);

        %PlotScale=0.1*max(abs(R(2*Nnodes+1:end)))*min([max(x)-min(x) max(y)-min(y)])/CtrlVar.PlotXYscale;
        PlotCircles(x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,R(2*Nnodes+1:end)/PlotScale,'r')
    end
    PlotGroundingLines(CtrlVar,MUA,F.GF) ;
    PlotCalvingFronts(CtrlVar,MUA,F);
    axis equal tight

    Ru=R(1:MUA.Nnodes);
    ru2=Ru'*Ru; 
    %FigRu=FindOrCreateFigure("Ru"); clf(FigRu) ;
    UaPlots(CtrlVar,MUA,F,Ru,GetRidOfValuesDownStreamOfCalvingFronts=false,FigureTitle="Ru residuals")  ;
    ModifyColormap ;
    title(sprintf("$R_u$ residuals. Iteration=%i. $r_u^2$=%g",iteration,ru2),Interpreter="latex")

    Rv=R(MUA.Nnodes+1:2*MUA.Nnodes);
    rv2=Rv'*Rv ;
    %FigRv=FindOrCreateFigure("Rv"); clf(FigRv) ;
    UaPlots(CtrlVar,MUA,F,Rv,GetRidOfValuesDownStreamOfCalvingFronts=false,FigureTitle="Rv residuals")  ;
    ModifyColormap ;
    title(sprintf("$R_v$ residuals. Iteration=%i. $r_v^2$=%g",iteration,rv2),Interpreter="latex")

    if contains(msg,'h')
        Rh=R(2*MUA.Nnodes+1:3*MUA.Nnodes);
        rh2=Rh'*Rh; 
        %FigRh=FindOrCreateFigure("Rh"); clf(FigRh) ;
        UaPlots(CtrlVar,MUA,F,Rh,GetRidOfValuesDownStreamOfCalvingFronts=false,FigureTitle="Rh residuals")  ;
        ModifyColormap ;
        title(sprintf("$R_h$ residuals. Iteration=%i. $r_h^2$=%g",iteration,rh2),Interpreter="latex")
    end


else % h residuals only


    hold on
    title('Nodal Force residuals (R+L^T \lambda)')

    PlotScale=0.1*max(abs(R))*min([max(x)-min(x) max(y)-min(y)])/CtrlVar.PlotXYscale;
    xmax=max(x) ; xmin=min(x) ; ymax=max(y) ;ymin=min(y);

    [R,I]=sort(abs(R),'descend') ; x=x(I) ; y=y(I);
    N=min([1000,numel(x)]);
    R=R(1:N) ; x=x(1:N) ; y=y(1:N);
    PlotCircles(x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,R/PlotScale,'r')
    axis([xmin xmax ymin ymax]/CtrlVar.PlotXYscale) ; axis equal tight

end
%%

drawnow
end

