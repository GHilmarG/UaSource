




function PlotResultsFromThicknessInversion(CtrlVar,MUA,F,BCs,Priors,Meas,hest,htrue)

narginchk(8,8)




UaPlots(CtrlVar,MUA,F,hest,FigureTitle="hest",CreateNewFigure=true) ;
title("h estimated")
hold on  ; PlotGroundingLines(); 



if ~isempty(Priors.h)
    UaPlots(CtrlVar,MUA,F,Priors.h,FigureTitle="hprior",CreateNewFigure=true) ;  title("h prior")
    hold on  ; PlotGroundingLines(); 
end

if ~isempty(Meas.h)
    
    UaPlots(CtrlVar,MUA,F,Meas.h,FigureTitle="hmeas",CreateNewFigure=true) ;  
    title("$h$ measured",Interpreter="latex")
    hold on  ; PlotGroundingLines(); 

    
    UaPlots(CtrlVar,MUA,F,hest-Meas.h,FigureTitle="hest-hmeas",CreateNewFigure=true) ;  
    title("$h$: estimated - measured ice thickness",Interpreter="latex")
    hold on  ; PlotGroundingLines(); 



end






if ~anynan(htrue) && ~isempty(htrue)

    UaPlots(CtrlVar,MUA,F,htrue,FigureTitle="h true") ; title("$h$ true",Interpreter="latex")
    axis tight

    FigCompare=FindOrCreateFigure("compare") ; clf(FigCompare)
    yyaxis left
    plot(htrue,hest,'.')
    hold on
    plot(xlim,xlim);
    ylim(xlim)
    ylabel("h estimated")
    yyaxis right
    plot(htrue,hest-htrue,'.')
    ylabel("hest-htrue")

    lm=fitlm(hest,htrue);
    title(sprintf("R2=%f",lm.Rsquared.Ordinary));
    xlabel("h true")  ;

    UaPlots(CtrlVar,MUA,F,hest-htrue,FigureTitle="hest-htrue",CreateNewFigure=true) ;  
    title("hest-htrue")

end

end