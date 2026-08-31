       
function [ForceFig,WorkFig]=PlotCostFunctionsVersusGamma(CtrlVar,RunInfo,gamma,r,iteration,IdString,...
        gammaTestVector,rForceTestvector,rWorkTestvector,rD2Testvector,...
        slopeForce,slopeWork,slopeD2,rForce,rWork,D2)
    
    [temp,I0]=min(abs(gammaTestVector)) ;
    Lower=min(gammaTestVector);
    Upper=max(gammaTestVector);
    
    ForceFig=FindOrCreateFigure(IdString+" Force and Work Residuals"+num2str(iteration));
    clf(ForceFig)
    hold off
    % slopeForce=-2*rForce0;
    yyaxis left
    plot(gammaTestVector,rForceTestvector,'o-',DisplayName="test vector force") ;
    hold on ;
    plot([gammaTestVector(I0) gammaTestVector(I0+1)],[rForceTestvector(I0) rForceTestvector(I0)+(gammaTestVector(I0+1)-gammaTestVector(I0))*slopeForce],'b-',LineWidth=2,...
        DisplayName="slope force")
    ylabel('Force Residuals')
    
    if CtrlVar.MinimisationQuantity=="Force Residuals"


        if size(RunInfo.BackTrack.Infovector,2)>1

            plot(RunInfo.BackTrack.Infovector(:,1),RunInfo.BackTrack.Infovector(:,2),"ok",MarkerFaceColor="k",DisplayName="Backtrack info vector force");
            plot(gamma,r,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r',MarkerSize=10,DisplayName="$(\gamma,r)$ Force Residuals")

        end

    end
    
    yyaxis right
    % slopeWork=-2*rWork0;
    plot(gammaTestVector,rWorkTestvector,'o-',DisplayName="test vector work") ; 
    hold on ;
    plot([gammaTestVector(I0) gammaTestVector(I0+1)],[rWorkTestvector(I0) rWorkTestvector(I0)+(gammaTestVector(I0+1)-gammaTestVector(I0))*slopeWork],'r-',...
        LineWidth=2,DisplayName="slope")
    ylabel('Work Residuals')
    
    if CtrlVar.MinimisationQuantity=="Work Residuals"
        plot(gamma,r,'Marker','h','MarkerEdgeColor','k','MarkerFaceColor','g',MarkerSize=2,DisplayName="$(\gamma,r)$ Work Residuals")
        
    end
    title(sprintf('%s iteration %-i,  iarm=%-i, $\\gamma$=%-6g , r=%-6g',IdString,iteration,RunInfo.BackTrack.iarm,gamma,r),'interpreter','latex') ;
    xlabel(' \gamma ') ;
    xlim([Lower Upper])
    lg=legend(Interpreter="latex");

    WorkFig=FindOrCreateFigure(IdString+" rWork and D2"+num2str(iteration));
    clf(WorkFig)
    hold off
    % slopeD2=-D20;
    yyaxis left
    plot(gammaTestVector,rD2Testvector,'o-',DisplayName="test vector D2") ; hold on ;
    plot([gammaTestVector(I0) gammaTestVector(I0+1)],[rD2Testvector(I0) rD2Testvector(I0)+(gammaTestVector(I0+1)-gammaTestVector(I0))*slopeD2],'b-',...
        LineWidth=2,DisplayName="slope D2")
    ylabel('$D^2$','interpreter','latex')
    % fig.CurrentAxes.XAxisLocation='origin';
    
    
    
    yyaxis right
    % slopeWork=-2*rWork0;
    plot(gammaTestVector,rWorkTestvector,'o-',DisplayName="Test vector work") ; 
    hold on ;
    
    plot([gammaTestVector(I0) gammaTestVector(I0+1)],[rWorkTestvector(I0) rWorkTestvector(I0)+(gammaTestVector(I0+1)-gammaTestVector(I0))*slopeWork],'r-',LineWidth=2,...
        DisplayName="Slope work") ; 
    ylabel('Work Residuals')
    
    if CtrlVar.MinimisationQuantity=="Work Residuals"
        
        plot(RunInfo.BackTrack.Infovector(:,1),RunInfo.BackTrack.Infovector(:,2),"ok",MarkerFaceColor="k",DisplayName="Backtrack Info Vector work"); 
        plot(gamma,r,'Marker','*','MarkerEdgeColor','r','MarkerFaceColor','r',DisplayName="$(\gamma,r)$ Work Residuals")
    end
    

    
    
    title(sprintf('%s iteration %-i,  iarm=%-i, $\\gamma$=%-6g, $r_W$=%-6.3g, $r_F$=%-6.3g',IdString,iteration,RunInfo.BackTrack.iarm,gamma,rWork,rForce),'interpreter','latex') ;
    xlabel(' \gamma ') ;
    WorkFig.CurrentAxes.XAxisLocation='origin'; grid on ;
    xlim([Lower Upper])
    hold off
    lg=legend(Interpreter="latex");
    
end