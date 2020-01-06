function [RunInfo,varargout]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin)
    
    %%
    % varargout=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,varargin)
    %
    % example:
    % [h2,s2,rho2]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,h1,s1,rho1)
    % interpolates h,s and rho, from FE Mesh MUA1 onto (x2,y2).
    
    %%
    
   
    nVarargsIn = length(varargin);
    nVar=nVarargsIn;
    varargout=cell(nVar,1);
    
    xOld=MUAold.coordinates(:,1);
    yOld=MUAold.coordinates(:,2);
    xNew=MUAnew.coordinates(:,1);
    yNew=MUAnew.coordinates(:,2);
    
    
    if isempty(xNew)
        for iOut=1:nVar
            varargout{iOut}=NaN;
        end
        return
    end
    
    
    
    % The underlying triangulation is only done once, and when interpolating different fields
    % only the values are changed.
    
    %tri=TriFE(MUA1.connectivity);  % Only corner points
    % tri=CreateFEmeshTriRep(MUA1.connectivity,MUA1.coordinates);
    
    %% First check if some (x2,y2) locations are identical to (x1,y1) locations
    % and where this is the case, just copy values across
    % TR = delaunayTriangulation(x1,y1);
    
    Tarea=TriAreaFE(MUAold.coordinates,MUAold.connectivity);
    tol=1e-5*sqrt(2*min(Tarea)) ;  % I found that tol=1000*eps is not enought...
    
    
    [ID,d] = nearestNeighbor(MUAold.TR,[xNew yNew]);  % This works for all element types! (3, 6 and 10)
    % The reason this works for all element types despite TR.ConnectivityList
    % only containing the corner nodes is because all coordinates are included
    % TR.Points.
    % For example this gives correct answer:
    %  T=triangulation([1 3 5],[0 0 ; 0.5 0 ; 1 0 ; 0.5 0.5 ; 0 1 ; 0 0.5]) ;
    %  [ID,d] = nearestNeighbor(T,[0.5 0]);
    
    
    IdenticalNodes=d<tol ;
    
    nNewNodes=numel(xNew);
    nOldNodes=numel(xOld);
    nIdenticalNodes=numel(find(IdenticalNodes));
    NotIdendicalNodes=find(~IdenticalNodes);
    nNotSame=numel(NotIdendicalNodes);
    
    RunInfo.Mapping.nNewNodes=nNewNodes;
    RunInfo.Mapping.nOldNodes=nOldNodes;
    RunInfo.Mapping.nIdenticalNodes=nIdenticalNodes;
    RunInfo.Mapping.nNotIdenticalNodes=nNotSame;
    
    if nNewNodes == nIdenticalNodes
        % all new nodes are identical to old ones
        RunInfo.Mapping.nNotIdenticalNodesOutside=0;
        RunInfo.Mapping.nNotIdenticalInside=0;
    end
    
    for iVar=1:nVar
        varargout{iVar}(IdenticalNodes)=varargin{iVar}(ID(IdenticalNodes));
    end
    
    
    
    %% Now check if any locations are remaining
    %
    
    if (nNewNodes-nIdenticalNodes)>0
        
        % Determine which of the new nodes are inside and which outside of the old
        % FE mesh over which a solution has already been calculated.
        %
        
        % For the remaining (x2,y2) locations:
        % 1) if outside values have been defined, use those values
        % 2) elsewere use scattered interpolant
        %
        
        % Are any of the remaining locations within the old mesh?
        
        ID = pointLocation(MUAold.TR,[xNew(~IdenticalNodes) yNew(~IdenticalNodes)]) ;  % ID is NaN if outside all triangles,
        
        NodesOutside=NotIdendicalNodes(isnan(ID)) ;
        NodesInsideAndNotSame=NotIdendicalNodes(~isnan(ID)) ;
        
        % Jan 2019:  poinLocation actually fails in some cases!
        % If a point is on an edge within a triangulation ID and barycentric
        % coordinates can be NaNs even if matlab documentation states that one of
        % the adjacent triangles will be selected.
        %
        % I guess one way of dealing with this is to add a small vector to the xy
        % locations in two orthonormal directions. This should work because I
        % already have delt with all idential nodes, so I now only have a potential
        % problem with xy locations along the edges of the triangulatino. If I
        % shift  in two orhtonormal directions I'm garanteed not to be moving along
        % the edge itself.
        %
        if ~isempty(NodesOutside)
            shift=tol;
            IDTest = pointLocation(MUAold.TR,[xNew(NodesOutside)+shift yNew(NodesOutside)+2*shift]) ;
            
            NewInsideAndNotSameNodes=NodesOutside(~isnan(IDTest));
            %  Add those nodes to the right set
            if ~isempty(NewInsideAndNotSameNodes)
                NodesInsideAndNotSame=[NodesInsideAndNotSame;NewInsideAndNotSameNodes];
                NodesOutside=setdiff(NodesOutside,NewInsideAndNotSameNodes);
            end
            
            IDTest = pointLocation(MUAold.TR,[xNew(NodesOutside)-2*shift yNew(NodesOutside)+shift]) ;
            NewInsideAndNotSameNodes=NodesOutside(~isnan(IDTest));
            %  Add those nodes to the right set
            if ~isempty(NewInsideAndNotSameNodes)
                NodesInsideAndNotSame=[NodesInsideAndNotSame;NewInsideAndNotSameNodes];
                NodesOutside=setdiff(NodesOutside,NewInsideAndNotSameNodes);
            end
            
            % And finally, are any of the outside nodes actually on the mesh boundary?
            NodesOnBoundary = DistanceToLineSegment([xNew(NodesOutside) yNew(NodesOutside)],[MUAold.Boundary.x MUAold.Boundary.y],[],1000*eps);
            
            %  Add any ouside nodes on bounday to the set of inside nodes
            if ~isempty(NodesOnBoundary)
                NodesInsideAndNotSame=[NodesInsideAndNotSame;NodesOutside(NodesOnBoundary)];
                NodesOutside(NodesOnBoundary)=[];
            end
        end
        
        
        RunInfo.Mapping.nNotIdenticalNodesOutside=numel(NodesOutside);
        RunInfo.Mapping.nNotIdenticalNodesInside=numel(NodesInsideAndNotSame);
        
        
        if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots
            FindOrCreateFigure("-Old and new nodes-");
            tt=axis;
            hold off
            p0=PlotMuaMesh(CtrlVar,MUAnew,[],'b');
            hold on
            
            p1=PlotMuaMesh(CtrlVar,MUAold,[],'k');
            if ~isequal(tt,[0 1 0 1])
                axis(tt)
            end
            
            p2=plot(xNew(NodesOutside)/CtrlVar.PlotXYscale,yNew(NodesOutside)/CtrlVar.PlotXYscale,'om');
            p3=plot(xNew(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,yNew(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,'or');
            p4=plot(xNew(IdenticalNodes)/CtrlVar.PlotXYscale,yNew(IdenticalNodes)/CtrlVar.PlotXYscale,'*g');
            
            if ~isempty(p2) && ~isempty(p3)
                legend([p0 p1 p2 p3 p4],'New Mesh','Old Mesh','New and outside','New but inside','Save','Location','northeastoutside')
            elseif ~isempty(p3)
                legend([p0 p1 p3 p4],'New Mesh','Old Mesh','New but inside','Same','Location','northeastoutside')
            elseif  ~isempty(p2)
                legend([p0 p1 p2 p4],'New Mesh','Old Mesh','New and outside','Same','Location','northeastoutside')
            end
            
            
            axis equal
            hold off
            
            
            
            
            % fprintf('#Outside=%i \t   #Inside and not same=%i \n',numel(NodesOutside),numel(NodesInsideAndNotSame))
            FigTitle=sprintf('             #Nodes in new mesh=%i \t #Nodes in old mesh=%i \t \n #Same Nodes=%i \t  #~Same Nodes=%i \t #Nodes inside and new=%i \t #Outside nodes=%i ',...
                nNewNodes,nOldNodes,nIdenticalNodes,nNotSame,numel(NodesInsideAndNotSame),numel(NodesOutside)) ;
            title(FigTitle)
        end
        
        F = scatteredInterpolant();
        F.Points=MUAold.coordinates;
        
        F.Method='natural';
        F.ExtrapolationMethod='nearest';
        
        % If OutsideValues have been defined for the variable, then use these for
        % the outside points. Otherwise use scattered interpolant.
        % For all the remaingin new nodes within the old mesh use scattered
        % interpolant
        
        for iVar=1:nVar
            
            F.Values=double(varargin{iVar});
            if ~isempty(OutsideValues) && ~isnan(OutsideValues(iVar))
                varargout{iVar}(NodesOutside)=OutsideValues(iVar);
            else
                varargout{iVar}(NodesOutside)=F(xNew(NodesOutside),yNew(NodesOutside));
            end
            varargout{iVar}(NodesInsideAndNotSame)=F(xNew(NodesInsideAndNotSame),yNew(NodesInsideAndNotSame));
        end
        
        
    end
    %% Now check that all returned variables are column vectors
    
    for iVar=1:nVar
        if ~isempty(varargout{iVar}) && ~iscolumn(varargout{iVar})
            varargout{iVar}=varargout{iVar}';
        end
    end
    
    
end