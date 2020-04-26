function [UserVar,LSF,BCsLevelSet,CalvingRate]=DefineCalving(UserVar,CtrlVar,MUA,BCs,BCsLevelSet,F)
    
    %%
    
    LSF=F.LSF ;
    
    if CtrlVar.CurrentRunStepNumber< 2
        
%       x=MUA.coordinates(:,1);
%       y=MUA.coordinates(:,2);
%       xMax=max(x) ; xMin=min(x) ;  xL=1e6*(xMax-xMin);
%       yMax=max(y) ; yMin=min(y) ;  yL=xL; 
%
% Various options to define the initial calving front:
%         % Option 1
%         CalvingFrontNodes=setdiff(MUA.Boundary.Nodes,union(BCs.ubFixedNode,BCs.vbFixedNode)) ; % CalvingFrontNodes
%         xc=MUA.coordinates(CalvingFrontNodes,1) ;
%         yc=MUA.coordinates(CalvingFrontNodes,2) ;
%         [yc,iSort]=sort(yc,'descend')  ; xc=xc(iSort) ;
%         CalvingFront=[xc(:{ yc(:)] ; 
%         
%         % Option 2
%         xc=600e3 ; 
%         CalvingFront=[xc  yMin ; xc yMax ] ; 
%         CalvingFrontClosure=[CalvingFront(end,1) yL ;   xMin-xL yL ; xMin-xL -yL ;  CalvingFront(1,1) -yL  ;CalvingFront(1,:)  ] ;
%         CalvingFront=[CalvingFront; CalvingFrontClosure ] ;
%         Npoints=10000 ; CalvingFront = interparc(Npoints,CalvingFront(:,1),CalvingFront(:,2),'linear'); % add some points
%         LSF=SignedDistance(MUA.coordinates,CalvingFront);

        % option 3
        xc=600e3; 
        LSF=xc-MUA.coordinates(:,1) ; 
        % figure ; plot(CalvingFront(:,1)/1000,CalvingFront(:,2)/1000,'-o'); axis equal
    end
    
    % I=MUA.coordinates(:,1) < 450e3 ; LSF(I)=10e4 ; % Possible re-initialisation 
    
    if CtrlVar.time<1
        CalvingRate=zeros(MUA.Nnodes,1); % always define the calving rate
    else
        CalvingRate=10*F.h;  
    end
    
    
    Delta = DiracDelta(1/UserVar.CalvingDeltaWidth,LSF,0) ; 
    Delta=Delta/max(Delta) ; 
    
    CalvingRate=CalvingRate.*Delta;
    
    %%
    % fig=FindOrCreateFigure('LSF in Define Calving') ; 
    % PlotMeshScalarVariable(CtrlVar,MUA,LSF) ; 
    
    %%
end
