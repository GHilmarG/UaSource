function  [UserVar,x,y,EleSizeDesired,EleSizeCurrent,ElementsToBeRefined, ElementsToBeCoarsened,NodalErrorIndicators]=...
    DesiredEleSizes(UserVar,CtrlVar,MUA,F,l,GF,Ruv,Lubvb)

%
% Estimates optimal element sizes based on number of explicit error estimators
% The estimates are given at the nodal points of the current FE mesh
%
% Scales element sizes to fit within the range of CtrlVar.MeshSizeMin to CtrlVar.MeshSizeMax
%

ubvbLambda=l.ubvb;

hf=(F.S-F.B)*F.rhow./F.rho ;

% Calculate current element sizes
EleArea=TriAreaFE(MUA.coordinates,MUA.connectivity);
M= Ele2Nodes(MUA.connectivity,MUA.Nnodes);
EleSizeCurrent=sqrt(M*EleArea);
% figure ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,EleSize0);
% title('Nodal areas')

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2); 

EleSizeDesired=zeros(numel(x),1)+CtrlVar.MeshSizeMax ;
NodalErrorIndicators=[];

for I=1:numel(CtrlVar.RefineCriteria)
    fprintf(CtrlVar.fidlog,' remeshing criterion is : %s \n ',CtrlVar.RefineCriteria{I});
    ErrorIndicatorUsefull=1;
    
    %% calculations specific to error criterion start
    switch CtrlVar.RefineCriteria{I}
        
        case 'effective strain rates'
            
            [~,~,~,~,~,~,~,e]=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n);
            NodalErrorIndicator=ProjectFintOntoNodes(MUA,e);
            NodalErrorIndicator(NodalErrorIndicator<0)=0;
            
            for II=1:CtrlVar.NumberOfSmoothingErrorIndicatorIterations
                EleErrorIndicator=Nodes2EleMean(MUA.connectivity,NodalErrorIndicator);  NodalErrorIndicator=M*EleErrorIndicator;
            end
            
            if ~isnan(CtrlVar.RefineCriteriaFlotationLimit(I))
                ind=abs(F.h-hf)>CtrlVar.RefineCriteriaFlotationLimit(I);
                NodalErrorIndicator(ind)=0;
            end
            
        % Using effective strain rates (e) for error esimation can tricky,
        % even if e is zero the errors there are not zero...
        % It is unlikely that this will happen in a real situation,
        % but sometimes e is very small for a limited number of elements.
        % When scaled up this gives distored estimates
        % I try to get around this problem by discarding the lowest 20% of estimates
        test=sort(NodalErrorIndicator);
        iThresh=round(numel(test)*0.2);
        if iThresh>1
            Thresh=test(iThresh);
            NodalErrorIndicator(NodalErrorIndicator<Thresh)=Thresh;
        end
            
            
            NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
            
            if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                figure(1710) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar); 
                               title('Relative error based on effective strain rates') ;
            end
            
            
            NodalErrorIndicators.StrainRates=NodalErrorIndicator;
 
 
        case '|dhdt|'
            
            NodalErrorIndicator=abs(F.dhdt);
            if all(NodalErrorIndicator<1e-5)  % this could for example happen at a start of a run where dhdt has not been calculated
                NodalErrorIndicator=NodalErrorIndicator*0 ;
                ErrorIndicatorUsefull=0;
            else
                % smooth this over a few elements
                
                for II=1:CtrlVar.NumberOfSmoothingErrorIndicatorIterations
                    EleErrorIndicator=Nodes2EleMean(MUA.connectivity,NodalErrorIndicator);  NodalErrorIndicator=M*EleErrorIndicator;
                end
                
                
                if ~isnan(CtrlVar.RefineCriteriaFlotationLimit(I))
                    ind=abs(F.h-hf)>CtrlVar.RefineCriteriaFlotationLimit(I);
                    NodalErrorIndicator(ind)=0;
                end
                
                
                NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
                NodalErrorIndicators.RateOfThicknessChange=NodalErrorIndicator;
                
                if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                    figure(1720) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);
                    title('Relative error based on |dhdt| ')
                end
            end
            
        case 'residuals'
            if ~isempty(Lubvb)
                Ruv=Ruv+Lubvb'*ubvbLambda;
            end
            
            NodalErrorIndicator=sqrt((Ruv(1:MUA.Nnodes)).^2+(Ruv(MUA.Nnodes+1:2*MUA.Nnodes)).^2);
            
            if all(NodalErrorIndicator<1e-5)  % this could for example happen at a start of a run where dhdt has not been calculated
                NodalErrorIndicator=NodalErrorIndicator*0 ;
                ErrorIndicatorUsefull=0;
            else
                % smooth this over a few elements
                
                for II=1:CtrlVar.NumberOfSmoothingErrorIndicatorIterations
                    EleErrorIndicator=Nodes2EleMean(MUA.connectivity,NodalErrorIndicator);  NodalErrorIndicator=M*EleErrorIndicator;
                end
                
                
                if ~isnan(CtrlVar.RefineCriteriaFlotationLimit(I))
                    ind=abs(F.h-hf)>CtrlVar.RefineCriteriaFlotationLimit(I);
                    NodalErrorIndicator(ind)=0;
                end
                
                
                NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
                
                NodalErrorIndicators.Residuals=NodalErrorIndicator;
                
                if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                    figure(1725) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);
                    title('Relative error based on residuals ')
                end
            end
            
        case '||grad(dhdt)||'
            
            if all(abs(dhdt)<1e-5)
                ErrorIndicatorUsefull=0;
                fprintf(CtrlVar.fidlog,' WARNING: dh/dt too small to be useful as an explicit error indicator \n ');
            else
                [dfdx,dfdy]=calcFEderivativesMUA(F.dhdt,MUA,CtrlVar);
                %dfdx=sum(dfdx,2) ; dfdy=sum(dfdy,2) ;
                EleErrorIndicator=sqrt(dfdx.*dfdx+dfdy.*dfdy);
                NodalErrorIndicator=ProjectFintOntoNodes(MUA,EleErrorIndicator);
                NodalErrorIndicator(NodalErrorIndicator<0)=0;
                for II=1:CtrlVar.NumberOfSmoothingErrorIndicatorIterations
                    EleErrorIndicator=Nodes2EleMean(MUA.connectivity,NodalErrorIndicator);  NodalErrorIndicator=M*EleErrorIndicator;
                end
                
                
                if ~isnan(CtrlVar.RefineCriteriaFlotationLimit(I))
                    ind=abs(F.h-hf)>CtrlVar.RefineCriteriaFlotationLimit(I);
                    NodalErrorIndicator(ind)=0;
                end
                
                NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
                NodalErrorIndicators.dhdtGradient=NodalErrorIndicator;
                
                if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                    figure(1730) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);
                    title('Relative error based on ||grad(dhdt)|| ')
                end
                
            end
            
        case 'dhdt curvature'
            
            if all(abs(dhdt)<1e-5)
                ErrorIndicatorUsefull=0;
                fprintf(CtrlVar.fidlog,' WARNING: dh/dt too small to be usefull as an explicit error indicator \n ');
            else
                
                [dfdx,dfdy]=calcFEderivativesMUA(F.dhdt,MUA,CtrlVar);
                %dfdx=sum(dfdx,2) ; dfdy=sum(dfdy,2) ;
                D=sqrt(1+dfdx.^2+dfdy.^2);
                fx=dfdx./D ; fy=dfdy./D;
                % project back onto nodes
                [fx,fy]=ProjectFintOntoNodes(MUA,fx,fy);
                [dfdxx,~]=calcFEderivativesMUA(fx,MUA,CtrlVar);
                [~,dfdyy]=calcFEderivativesMUA(fy,MUA,CtrlVar);
                %dfdxx=sum(dfdxx,2) ; dfdyy=sum(dfdyy,2) ;
                
                EleErrorIndicator=(abs(dfdxx)+abs(dfdyy))/2;
                                
                NodalErrorIndicator=ProjectFintOntoNodes(MUA,EleErrorIndicator);
                NodalErrorIndicator(NodalErrorIndicator<0)=0;
                %NodalErrorIndicator=M*EleErrorIndicator;
                %
                %                     R=1./NodalErrorIndicator;
                %                     fprintf('Inverse curvature (R) based on dhdt is: max(R)=%-g \t min(R)=%-g \t median(R)=%-g \n ',max(R),min(R),median(R))
                %                     fprintf('Min of inverse curvature suggests an element size scale of considerably less than %-g. Min prescribed ele size is %-g \n',min(R),CtrlVar.MeshSizeMin);
                %                     fprintf('The ratio: min(R)/MeshSizeMin=%-g \n ',min(R)/CtrlVar.MeshSizeMin);
                
                for II=1:CtrlVar.NumberOfSmoothingErrorIndicatorIterations
                    EleErrorIndicator=Nodes2EleMean(MUA.connectivity,NodalErrorIndicator);  NodalErrorIndicator=M*EleErrorIndicator;
                end
                
                
                if ~isnan(CtrlVar.RefineCriteriaFlotationLimit(I))
                    ind=abs(F.h-hf)>CtrlVar.RefineCriteriaFlotationLimit(I);
                    NodalErrorIndicator(ind)=0;
                end
                
                NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
                
                NodalErrorIndicators.dhdtCurvature=NodalErrorIndicator;
                
                if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                    figure(1740) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);
                    title('Relative error based on dhdt curvature ')
                end
                
            end
        case 'thickness gradient'
            
            if (max(F.h)-min(F.h))< 10
                ErrorIndicatorUsefull=0;
                fprintf(CtrlVar.fidlog,' WARNING: Thickness variation too small to be usefull as an explicit error indicator \n ');
            else
                [dfdx,dfdy]=calcFEderivativesMUA(F.h,MUA,CtrlVar);
                %dfdx=sum(dfdx,2) ; dfdy=sum(dfdy,2) ;
                EleErrorIndicator=sqrt(dfdx.*dfdx+dfdy.*dfdy);
                NodalErrorIndicator=ProjectFintOntoNodes(MUA,EleErrorIndicator);
                NodalErrorIndicator(NodalErrorIndicator<0)=0;               
                for II=1:CtrlVar.NumberOfSmoothingErrorIndicatorIterations
                    EleErrorIndicator=Nodes2EleMean(MUA.connectivity,NodalErrorIndicator);  
                    NodalErrorIndicator=M*EleErrorIndicator;
                end
                
                
                if ~isnan(CtrlVar.RefineCriteriaFlotationLimit(I))
                    ind=abs(F.h-hf)>CtrlVar.RefineCriteriaFlotationLimit(I);
                    NodalErrorIndicator(ind)=0;
                end
                
                NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
                NodalErrorIndicators.ThicknessGradient=NodalErrorIndicator;
                if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                    figure(1750) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);
                    title('Relative error based on norm of thickness gradient ')
                end
                
            end
            
        case 'thickness curvature'
            
            if (max(F.h)-min(F.h))< 10
                ErrorIndicatorUsefull=0;
                fprintf(CtrlVar.fidlog,' WARNING: Thickness variation too small to be usefull as an explicit error indicator \n ');
            else
                
                [dfdx,dfdy,]=calcFEderivativesMUA(F.h,MUA,CtrlVar);
               % dfdx=sum(dfdx,2) ; dfdy=sum(dfdy,2) ;
                D=sqrt(1+dfdx.^2+dfdy.^2);
                dfdx=dfdx./D ; dfdy=dfdy./D;
                % project back onto nodes
                [dfdx,dfdy]=ProjectFintOntoNodes(MUA,dfdx,dfdy);
                
                [dfdxx,~]=calcFEderivativesMUA(dfdx,MUA,CtrlVar);
                [~,dfdyy]=calcFEderivativesMUA(dfdy,MUA,CtrlVar);
               % dfdxx=sum(dfdxx,2) ; dfdyy=sum(dfdyy,2) ;
                EleErrorIndicator=(abs(dfdxx)+abs(dfdyy))/2;
                %EleErrorIndicator=sqrt(dfdxx.*dfdxx+dfdxy.*dfdxy+dfdyx.*dfdyx+dfdyy.*dfdyy);
                
                
                NodalErrorIndicator=ProjectFintOntoNodes(MUA,EleErrorIndicator);
                NodalErrorIndicator(NodalErrorIndicator<0)=0;
                R=1./NodalErrorIndicator;
                fprintf('Inverse curvature (R) based on h is: max(R)=%-g \t min(R)=%-g \t median(R)=%-g \n ',max(R),min(R),median(R))
                fprintf('Min of inverse curvature suggests an element size scale of considerably less than %-g. Min prescribed ele size is %-g \n',min(R),CtrlVar.MeshSizeMin);
                fprintf('The ratio: min(R)/MeshSizeMin=%-g \n ',min(R)/CtrlVar.MeshSizeMin);
                
                
                for II=1:CtrlVar.NumberOfSmoothingErrorIndicatorIterations
                    EleErrorIndicator=Nodes2EleMean(MUA.connectivity,NodalErrorIndicator);  NodalErrorIndicator=M*EleErrorIndicator;
                end
                
                if ~isnan(CtrlVar.RefineCriteriaFlotationLimit(I))
                    ind=abs(F.h-hf)>CtrlVar.RefineCriteriaFlotationLimit(I);
                    NodalErrorIndicator(ind)=0;
                end
                
                NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
                NodalErrorIndicators.ThicknessCurvature=NodalErrorIndicator;
                if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                    figure(1760) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);
                    title('Relative error based on norm of thickness curvature ')
                end
                
            end
            
        case 'flotation'
            
            dgf = DiracDelta(1/CtrlVar.RefineDiracDeltaWidth,F.h-hf,CtrlVar.RefineDiracDeltaOffset);
            dgfmin=max(dgf)/1e5; dgf(dgf<dgfmin)=dgfmin;
            NodalErrorIndicator=dgf;
            
            NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
            
            NodalErrorIndicators.Flotation=NodalErrorIndicator;
            
            if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                figure(1770) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);  title(' error estimate based on flotation criterion')
                figure(1780) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,F.h-hf,CtrlVar);                 title(' h-hf ') ;
                figure(1790) ; plot(F.h-hf,dgf,'.'); title('Delta-function error indicator as a function of h-hf')
            end
            
            
            
        case 'f factor'
            
            k=1/CtrlVar.RefineDiracDeltaWidth;
            NodalErrorIndicator=(1.-exp(-k*(F.h-hf)))./(1+exp(-k*(F.h-hf)));
            
            NodalErrorIndicator=NodalErrorIndicator/max(NodalErrorIndicator);
            
            NodalErrorIndicators.fFactor=NodalErrorIndicator;
            
            if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
                figure(1800) ; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,NodalErrorIndicator,CtrlVar);
                title('Relative error based on f factor ')
            end
            
        case ''
                
           ErrorIndicatorUsefull=0;
            
        otherwise
            
            fprintf(CtrlVar.fidlog,' remeshing criterion is : %s \n ',CtrlVar.RefineCriteria{I});
            error(' what case? ')
            
    end
    
    %% calculations specific to error criterion done
    
    if ~ErrorIndicatorUsefull
        EleSizeIndicator =  zeros(numel(x),1)+CtrlVar.MeshSizeMax;
        fprintf(CtrlVar.fidlog,' Error indicator too similar across the mesh to be usefull. \n ');
    else
        
        
        
        
        % large explicit error indicator suggest small elements.
        
 
        % There is a potential numerical issue with the ErrorIndicator being exceedingly small for some elements
        % so in order to avoid this I define a smallest error value as:
        SmallestValue=max(NodalErrorIndicator)*CtrlVar.MeshSizeMin/CtrlVar.MeshSizeMax;
        NodalErrorIndicator(NodalErrorIndicator<SmallestValue)=SmallestValue;
        
        % now the error indicator is mapped into a ele size indicator
        % here some theory is needed, in general error is expected to
        % decrease inversely with ele size to some power
        EleSizeIndicator=real(1./((NodalErrorIndicator+SmallestValue).^CtrlVar.hpower)) ;
        
        % Ele=Ele0*Error
        % max=1/(err0)^p  -> err0=1/max^p
        % Now scale the error indicator to the whole range desired ele sizes
        % (this could be improved for example by weighting it with area)
        
        % I now do a simple weighting of each mesh size criterion
        % if I change elements sizes over some upper-fraction of the whole user-defined
        % range CtrlVar.MeshSizeMin to CtrlVar.MeshSizeMax
        % The fraction is given by CtrlVar.RefineCriteriaWeights(I)
        % Example:
        % If CtrlVar.RefineCriteriaWeights(I)=1 I used the whole range
        % If CtrlVar.RefineCriteriaWeights(I)=0.5 element size will range from
        % CtrlVar.MeshSizeMax down to CtrlVar.MeshSizeMin+(CtrlVar.MeshSizeMax-CtrlVar.MeshSizeMin)*(1-RefineCriteriaWeight)
        % If CtrlVar.RefineCriteriaWeights(I)=0 the criterion is effectivly ignored
        %
        %
        MeshSizeMax=CtrlVar.MeshSizeMax;
        
        MeshSizeMin=CtrlVar.MeshSizeMin+(CtrlVar.MeshSizeMax-CtrlVar.MeshSizeMin)*(1-CtrlVar.RefineCriteriaWeights(I));
        
        
        EleSizeIndicator=MeshSizeMin+...
            (MeshSizeMax-MeshSizeMin)*(EleSizeIndicator-min(EleSizeIndicator))/(max(EleSizeIndicator)-min(EleSizeIndicator));
        
    end
    
    
    if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
        figure(1810) ; hold off ;
        PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,EleSizeIndicator,CtrlVar);
        title(sprintf(' desired ele sizes %s ',CtrlVar.RefineCriteria{I}))
    end
    
    % take the minimum of each error indicator as measure for ele size
    EleSizeDesired=min(EleSizeDesired,EleSizeIndicator);
    
end


if all(EleSizeDesired==CtrlVar.MeshSizeMax)
    fprintf(' After using relative error criteria, all desired ele sizes are equal to CtrlVar.MeshSizeMax=%g.\n',CtrlVar.MeshSizeMax) 
    fprintf(' This most likely happened because either no relative error criteria were specified, or none were applicable. \n')
    fprintf(' Will now set all ele sizes equal to CtrlVar.MeshSize=%g \n ',CtrlVar.MeshSize);
    
    EleSizeDesired=zeros(MUA.Nnodes,1)+CtrlVar.MeshSize;
end


% do not allow EleSize to change too much and take a weighted average of
% the previous and the new EleSize:
EleSizeDesired=0.95*EleSizeDesired+0.05*EleSizeCurrent;

% and also put strickt limits on change in EleSize:
EleSizeRatio=EleSizeDesired./EleSizeCurrent;

I=EleSizeRatio>CtrlVar.MaxRatioOfChangeInEleSizeDuringAdaptMeshing; 
EleSizeDesired(I)=CtrlVar.MaxRatioOfChangeInEleSizeDuringAdaptMeshing*EleSizeCurrent(I);

I=EleSizeRatio<CtrlVar.MinRatioOfChangeInEleSizeDuringAdaptMeshing ; 
EleSizeDesired(I)=CtrlVar.MinRatioOfChangeInEleSizeDuringAdaptMeshing*EleSizeCurrent(I);


if CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
    figure(1820) ;
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,EleSizeDesired,CtrlVar);
    colorbar ;  title(sprintf(' combined desired ele sizes ')) ;
end


%%
% Set elesizes around GL to a specific value


if isfield(CtrlVar,'MeshAdapt') && isfield(CtrlVar.MeshAdapt,'GLrange')
    
    fprintf('Remeshing based on distance of nodes from grounding line.\n')
    
    KdTree=[];
    CtrlVar.PlotGLs=0;
    CtrlVar.GLsubdivide=1;
    [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,GF);
    
    
    
    
    for I=1:size(CtrlVar.MeshAdapt.GLrange,1)
        
        ds=CtrlVar.MeshAdapt.GLrange(I,1);
        dh=CtrlVar.MeshAdapt.GLrange(I,2);
        if dh<CtrlVar.MeshSizeMin
            fprintf('---> Warning: CtrlVar.MeshAdapt.GLrange(%i,2)=%g<CtrlVar.MeshSizeMin=%g \n',I,dh,CtrlVar.MeshSizeMin)
            fprintf('              Setting CtrlVar.MeshAdapt.GLrange(%i,2)=%g \n',I,CtrlVar.MeshSizeMin)
            dh=CtrlVar.MeshSizeMin;
        end
        fprintf('Nodes within the distance of %g from the grounding line are given the target element size %g \n',ds,dh)
        
        ID=FindAllNodesWithinGivenRangeFromGroundingLine(CtrlVar,MUA,xGL,yGL,ds,KdTree);
        
        EleSizeIndicator(ID)=dh;
        EleSizeDesired=min(EleSizeDesired,EleSizeIndicator);
    end
    
    
end
%%



% No further user defined modifications to EleSizeDesired

%% now create a (logical) list of elements to be refined, i.e. split up into four sub-elements

switch CtrlVar.RefineCriteria{1}
    
    case 'residuals'
        
        EleResiduals=Nodes2EleMean(MUA.connectivity,NodalErrorIndicators.Residuals);
        EleResidualsSorted=sort(EleResiduals);
        ElementsToBeRefined=EleResiduals>EleResidualsSorted(end-CtrlVar.LocalAdaptMeshMaxNrOfElementsToBeRefined);
        ElementsToBeCoarsened=false(MUA.Nele,1);
    otherwise
        
        eRatio=EleSizeDesired./EleSizeCurrent;
        eRatio=Nodes2EleMean(MUA.connectivity,eRatio);
        
        % do not refine a greater number of elements than CtrlVar.MeshRefinementRatio*CurrentNumberOfElements
        % at any given refinement step
        
        test=sort(eRatio);
        ElementsToBeRefined=eRatio<=test(ceil(numel(eRatio)*CtrlVar.LocalAdaptMeshRatio)) & eRatio<1;
        
        % have to make sure that if an element has just been refined that it will not
        % then afterwards be a candidate for coarsening. If an element was refined, the
        % size decreased by about a factor of 2 so if the ratio was 1+eps it is now
        % 0.5+eps and I must set eRatio>2 at the very least, for coarsening
        
        ElementsToBeCoarsened=eRatio>=test(floor(numel(eRatio)*CtrlVar.LocalAdaptMeshRatio)) & eRatio>4;
end



%% Now finally a user modification to EleSizeDesired and ElementsToBeRefined



[UserVar,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=GetDesiredEleSize(UserVar,CtrlVar,MUA,F,GF,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators);

assert(numel(x)==numel(y) && numel(x)==numel(EleSizeDesired),' Number of elements in x, y, and EleSize must be the same \n')


if   CtrlVar.doplots==1 && CtrlVar.doAdaptMeshPlots==1 && CtrlVar.InfoLevelAdaptiveMeshing>=10
    figure(1830)
    subplot(1,2,1,'replace')
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,EleSizeDesired,CtrlVar);
    title(' final desired ele sizes after user modification ');

    subplot(1,2,2,'replace')
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,EleSizeCurrent,CtrlVar);
    title(' current ele sizes  ');
    hold off
end




if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
    
    %figure ; PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,double(ElementsToBeRefined))  ;  title(' Refine Elements')
    
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    figure ;
    PlotMuaMesh(CtrlVar,MUA,[],'k');
    hold on
    PlotMuaMesh(CtrlVar,MUA,ElementsToBeRefined,'b');
    PlotMuaMesh(CtrlVar,MUA,ElementsToBeCoarsened,'r');
    
    
    title('Elements to be refined/coarsened in blue/red')
    fprintf('  Number of elements to be refined: %i \n',numel(find(ElementsToBeRefined)))
    fprintf('Number of elements to be coarsened: %i \n',numel(find(ElementsToBeCoarsened)))
    
end




end






