function [UserVar,LSF,BCsLevelSet,CalvingRate]=DefineCalving(UserVar,CtrlVar,MUA,BCs,BCsLevelSet,F)
    
    %%  This is the only DefineCalving m-file at the moment, still testing
    
    LSF=F.LSF ;
    
    if CtrlVar.CurrentRunStepNumber< 2

        xc=93.1e3; 
        LSF=xc-MUA.coordinates(:,1) ;
        % figure ; plot(CalvingFront(:,1)/1000,CalvingFront(:,2)/1000,'-o'); axis equal
    end
    
    
    if CtrlVar.time<0.5
        CalvingRate=zeros(MUA.Nnodes,1); % always define the calving rate
    else

        q=-2; 
        % k=u(350km)/h(350km)^q ; 
        k=86992392.1760846; 
        CalvingRate=k*F.h.^q;
        
        CalvingRate=CalvingRate(:) ;
        CalvingRate=CalvingRate.*(1-F.GF.node) ;
  
    end
    
    
    
 
    
    
%     Delta = DiracDelta(1/UserVar.CalvingDeltaWidth,LSF,0) ; 
%     Delta=Delta/max(Delta) ; 
%     
%     CalvingRate=CalvingRate.*Delta;
%     
    %%
    % fig=FindOrCreateFigure('LSF in Define Calving') ; 
    % PlotMeshScalarVariable(CtrlVar,MUA,LSF) ; 
    
    %%
end
