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
    
    x1=MUAold.coordinates(:,1);
    y1=MUAold.coordinates(:,2);
    x2=MUAnew.coordinates(:,1);
    y2=MUAnew.coordinates(:,2);
    
    
    if isempty(x2)
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
    tol=1e-5*sqrt(2*min(Tarea)) ;
  
    
    [ID,d] = nearestNeighbor(MUAold.TR,[x2 y2]);  % This works for all element types! (3, 6 and 10)
    same=d<tol ;
    
    % The reason this works for all elemetn types despite TR.ConnectivityList
    % only containing the corner nodes is because all coordinates are included
    % TR.Points.
    % For example this gives correct answer:
    %  T=triangulation([1 3 5],[0 0 ; 0.5 0 ; 1 0 ; 0.5 0.5 ; 0 1 ; 0 0.5]) ;
    %  [ID,d] = nearestNeighbor(T,[0.5 0]);
    %
    %
    
  
    
    
    nNew=numel(x2);
    nOld=numel(x1);
    nSame=numel(find(same));
    NodesNotSame=find(~same);
    nNotSame=numel(NodesNotSame);

    RunInfo.Mapping.nNewNodes=nNew;
    RunInfo.Mapping.nOldNodes=nOld;
    RunInfo.Mapping.nIdenticalNodes=nSame;
    RunInfo.Mapping.nNotIdenticalNodes=nNotSame;
    
    
    for iVar=1:nVar
        varargout{iVar}(same)=varargin{iVar}(ID(same));
    end
    
    
    
    %% Now check if any locations are remaining
    %
    
    if (nNew-nSame)>0
        
        % Determine which of the new nodes are inside and which outside of the old
        % FE mesh over which a solution has already been calculated. 
        %
        
        % For the remaining (x2,y2) locations:
        % 1) if outside values have been defined, use those values
        % 2) elsewere use scattered interpolant
        %
        
        % Are any of the remaining locations within the old mesh?
        
        ID = pointLocation(MUAold.TR,[x2(~same) y2(~same)]) ;  % ID is NaN if outside all triangles, 
        
        NodesOutside=NodesNotSame(isnan(ID)) ;
        NodesInsideAndNotSame=NodesNotSame(~isnan(ID)) ;
        
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
            IDTest = pointLocation(MUAold.TR,[x2(NodesOutside)+shift y2(NodesOutside)+2*shift]) ;
            
            NewInsideAndNotSameNodes=NodesOutside(~isnan(IDTest));
            %  Add those nodes to the right set
            if ~isempty(NewInsideAndNotSameNodes)
                NodesInsideAndNotSame=[NodesInsideAndNotSame;NewInsideAndNotSameNodes];
                NodesOutside=setdiff(NodesOutside,NewInsideAndNotSameNodes);
            end
            
            IDTest = pointLocation(MUAold.TR,[x2(NodesOutside)-2*shift y2(NodesOutside)+shift]) ;
            NewInsideAndNotSameNodes=NodesOutside(~isnan(IDTest));
            %  Add those nodes to the right set
            if ~isempty(NewInsideAndNotSameNodes)
                NodesInsideAndNotSame=[NodesInsideAndNotSame;NewInsideAndNotSameNodes];
                NodesOutside=setdiff(NodesOutside,NewInsideAndNotSameNodes);
            end
            
            % And finally, are any of the outside nodes actually on the mesh boundary?
            NodesOnBoundary = DistanceToLineSegment([x2(NodesOutside) y2(NodesOutside)],[MUAold.Boundary.x MUAold.Boundary.y],[],1000*eps);
            
            %  Add any ouside nodes on bounday to the set of inside nodes
            if ~isempty(NodesOnBoundary)
                NodesInsideAndNotSame=[NodesInsideAndNotSame;NodesOutside(NodesOnBoundary)];
                NodesOutside(NodesOnBoundary)=[];
            end
        end
 
     
        RunInfo.Mapping.nNotIdenticalNodesOutside=numel(NodesOutside);
        RunInfo.Mapping.nNotIdenticalInside=numel(NodesInsideAndNotSame);
        
        
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
            
            p2=plot(x2(NodesOutside)/CtrlVar.PlotXYscale,y2(NodesOutside)/CtrlVar.PlotXYscale,'om');
            p3=plot(x2(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,y2(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,'or');
            p4=plot(x2(same)/CtrlVar.PlotXYscale,y2(same)/CtrlVar.PlotXYscale,'*g');
            
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
                nNew,nOld,nSame,nNotSame,numel(NodesInsideAndNotSame),numel(NodesOutside)) ;
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
                varargout{iVar}(NodesOutside)=F(x2(NodesOutside),y2(NodesOutside));
            end
            varargout{iVar}(NodesInsideAndNotSame)=F(x2(NodesInsideAndNotSame),y2(NodesInsideAndNotSame));
        end
        
        
    end
    %% Now check that all returned variables are column vectors
    
    for iVar=1:nVar
        if ~isempty(varargout{iVar}) && ~iscolumn(varargout{iVar})
            varargout{iVar}=varargout{iVar}';
        end
    end
    
    
end