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
    
    isMapped=false(MUAnew.Nnodes,1); 
    
    
    % The underlying triangulation is only done once, and when interpolating different fields
    % only the values are changed.
    
    %tri=TriFE(MUA1.connectivity);  % Only corner points
    % tri=CreateFEmeshTriRep(MUA1.connectivity,MUA1.coordinates);
    
    %% First check if some (x2,y2) locations are identical to (x1,y1) locations
    % and where this is the case, just copy values across
    % TR = delaunayTriangulation(x1,y1);
    
    Tarea=TriAreaFE(MUAold.coordinates,MUAold.connectivity);
    tol=1e-5*sqrt(2*min(Tarea)) ;  % I found that tol=1000*eps is not enough...
    
    if isempty(MUAold.TR)
        MUAold.TR=CreateFEmeshTriRep(MUAold.connectivity,MUAold.coordinates);
    end

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
    nNotIdendicalNodes=numel(NotIdendicalNodes);
    
    RunInfo.Mapping.nNewNodes=nNewNodes;
    RunInfo.Mapping.nOldNodes=nOldNodes;
    RunInfo.Mapping.nIdenticalNodes=nIdenticalNodes;
    RunInfo.Mapping.nNotIdenticalNodes=nNotIdendicalNodes;
    
    if nNewNodes == nIdenticalNodes
        % all new nodes are identical to old ones
        RunInfo.Mapping.nNotIdenticalNodesOutside=0;
        RunInfo.Mapping.nNotIdenticalInside=0;
    end
    

    for iVar=1:nVar
        if isempty(varargin{iVar})
            varargout{iVar}=[];
        else
            varargout{iVar}(IdenticalNodes)=varargin{iVar}(ID(IdenticalNodes));
            isMapped(IdenticalNodes)=true;  % this could be outside the loop I think, and only needs to be done once as it is only dependedn on the values of IdenticalNodes
        end
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
            NodesOnBoundary = DistanceToLineSegment([xNew(NodesOutside) yNew(NodesOutside)],[MUAold.Boundary.x MUAold.Boundary.y],[],tol);
            
            %  Add any ouside nodes on bounday to the set of inside nodes
            if ~isempty(NodesOnBoundary)
                NodesInsideAndNotSame=[NodesInsideAndNotSame;NodesOutside(NodesOnBoundary)];
                NodesOutside(NodesOnBoundary)=[];
            end
        end
        
        
        RunInfo.Mapping.nNotIdenticalNodesOutside=numel(NodesOutside);
        RunInfo.Mapping.nNotIdenticalNodesInside=numel(NodesInsideAndNotSame);
        
        
        if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=5 
            fig=FindOrCreateFigure("-Old and new nodes-"); clf(fig) ; 
           
            tt=axis;
            hold off
            p0=PlotMuaMesh(CtrlVar,MUAnew,[],'b');
            hold on
            
            p1=PlotMuaMesh(CtrlVar,MUAold,[],'k');
            if ~isequal(tt,[0 1 0 1])
                axis(tt)
            end
            
            p2=plot(xNew(NodesOutside)/CtrlVar.PlotXYscale,yNew(NodesOutside)/CtrlVar.PlotXYscale,marker="h",color="m");
            p3=plot(xNew(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,yNew(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,'or');
            p4=plot(xNew(IdenticalNodes)/CtrlVar.PlotXYscale,yNew(IdenticalNodes)/CtrlVar.PlotXYscale,'*g');
            
            if ~isempty(p2) && ~isempty(p3)
                legend([p0 p1 p2 p3 p4],'New Mesh','Old Mesh','New and outside','New but inside','Identical','Location','northeastoutside')
            elseif ~isempty(p3)
                legend([p0 p1 p3 p4],'New Mesh','Old Mesh','New but inside','Identical','Location','northeastoutside')
            elseif  ~isempty(p2)
                legend([p0 p1 p2 p4],'New Mesh','Old Mesh','New and outside','Identical','Location','northeastoutside')
            end
            
            
            axis equal
            hold off
            
            
            
            
            % fprintf('#Outside=%i \t   #Inside and not same=%i \n',numel(NodesOutside),numel(NodesInsideAndNotSame))
            FigTitle=sprintf('             #Nodes in new mesh=%i \t #Nodes in old mesh=%i \t \n #Same Nodes=%i \t  #~Same Nodes=%i \t #Nodes inside and new=%i \t #Outside nodes=%i ',...
                nNewNodes,nOldNodes,nIdenticalNodes,nNotIdendicalNodes,numel(NodesInsideAndNotSame),numel(NodesOutside)) ;
            title(FigTitle)
        end
        

        %% Here for the first time I use scatteredInterpolant. This is only used for nodes that are not identical to those of the old 
        Finterpolant = scatteredInterpolant();
        Finterpolant.Points=MUAold.coordinates;
        
        Finterpolant.Method='natural';
        Finterpolant.ExtrapolationMethod='nearest';
        
        % If OutsideValues have been defined for the variable, then use these for
        % the outside points. Otherwise use scattered interpolant.
        % For all the remaining new nodes within the old mesh use scattered
        % interpolant
        
        for iVar=1:nVar
            
            if isempty(varargin{iVar})
                varargout{iVar}=[];
            else
                
                Finterpolant.Values=double(varargin{iVar});
                if ~isempty(OutsideValues) && ~isnan(OutsideValues(iVar))
                    varargout{iVar}(NodesOutside)=OutsideValues(iVar);
                    isMapped(NodesOutside)=true;
                else
                    varargout{iVar}(NodesOutside)=Finterpolant(xNew(NodesOutside),yNew(NodesOutside));
                    isMapped(NodesOutside)=true;
                end
                varargout{iVar}(NodesInsideAndNotSame)=Finterpolant(xNew(NodesInsideAndNotSame),yNew(NodesInsideAndNotSame));
                isMapped(NodesInsideAndNotSame)=true;
            end
        end
        
        
    end
    %% Now check that all returned variables are column vectors
    
    for iVar=1:nVar
        if ~isempty(varargout{iVar}) && ~iscolumn(varargout{iVar})
            varargout{iVar}=varargout{iVar}';
        end
    end
    

    if ~all(isMapped)

        find(~isMapped)
        error("MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant:NoAllMapped","not all nodes of the new mesh have been interpolated from the old one....!")

    end

    
end