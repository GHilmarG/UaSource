function [UserVar,LSF,BCsLevelSet,CalvingRate]=DefineCalving(UserVar,CtrlVar,MUA,BCs,BCsLevelSet,F)
    
    %%
    
    LSF=F.LSF ;
    
    if CtrlVar.CurrentRunStepNumber< 2
        
        x=MUA.coordinates(:,1);
        y=MUA.coordinates(:,2);
        
        xMax=max(x) ; xMin=min(x) ;  xL=xMax-xMin;
        yMax=max(y) ; yMin=min(y) ;  yL=yMax-yMin;
        
        % Various options to define the initial calving front:
%         % Option 1
%         CalvingFrontNodes=setdiff(MUA.Boundary.Nodes,union(BCs.ubFixedNode,BCs.vbFixedNode)) ; % CalvingFrontNodes
%         xc=MUA.coordinates(CalvingFrontNodes,1) ;
%         yc=MUA.coordinates(CalvingFrontNodes,2) ;
%         [yc,iSort]=sort(yc,'descend')  ; xc=xc(iSort) ;
%         CalvingFront=[xc(:{ yc(:)] ; 
%         
        % Option 2
        CalvingFront=[600e3  yMin ; 600e3 yMax ] ; 
            
   
        
        CalvingFrontClosure=[CalvingFront(end,1) yMax+yL ;   xMin-xL yMax+yL ; xMin-xL yMin-yL ;  CalvingFront(1,1) yMin-yL ] ;
        
        CalvingFront=[CalvingFront; CalvingFrontClosure ] ;
        
        Npoints=1000 ; CalvingFront = interparc(Npoints,CalvingFront(:,1),CalvingFront(:,2),'linear'); % add some points
        LSF=SignedDistance(MUA.coordinates,CalvingFront);
        
    end
    
    
    CalvingRate=-100e3 ; % always define the calving rate
    
    
    
    %%
    
    % figure ; PlotMeshScalarVariable(CtrlVar,MUA,LSF) ; colorbar
    
    %%
end