function varargout=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,MUA1,x2,y2,OutsideValues,varargin)
    
    %%
    % varargout=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,varargin)
    %
    % example:
    % [h2,s2,rho2]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,h1,s1,rho1)
    % interpolates h,s and rho, from FE Mesh MUA1 onto (x2,y2).
    
    %%
    
    nVarargsIn = length(varargin);
    
    %fprintf('MapNodalVariablesFromMesh1ToMesh2 : Total number of varargin inputs = %d\n',nVarargsIn);
    
    nVar=nVarargsIn;
    
    varargout=cell(nVar,1);
    
    if isempty(x2)
        for iOut=1:nVar
            varargout{iOut}=NaN;
        end
        return
    end
    
    x1=MUA1.coordinates(:,1);
    y1=MUA1.coordinates(:,2);
    x2=x2(:);
    y2=y2(:);
    
    
    
    % The underlying triangulation is only done once, and when interpolating different fields
    % only the values are changed.
    
    %tri=TriFE(MUA1.connectivity);  % Only corner points
    % tri=CreateFEmeshTriRep(MUA1.connectivity,MUA1.coordinates);
    
    %% First check if some (x2,y2) locations are identical to (x1,y1) locations
    % and where this is the case, just copy values across
    % TR = delaunayTriangulation(x1,y1);
    
    tol = eps*1000;
    [ID,d] = nearestNeighbor(MUA1.TR,[x2 y2]);  % only complete for 3-node elements
    same=d<tol ;
    
    %%  [ test 
    Tarea=TriAreaFE(MUA1.coordinates,MUA1.connectivity); Tlength=sqrt(2*min(Tarea))/100 ; 
    TestSame=d<Tlength ;
    
    if ~isequal(same,TestSame)
        save TestSave
        error(' testing TestSave ')
    end
    %% test ]
    
    
    nNew=numel(x2);
    nOld=numel(x1);
    nSame=numel(find(same));
    nNotSame=numel(find(~same));
    fprintf('#New=%i \t   #Old=%i \t  #Same=%i \t  #NotSame=%i \t #(New-Same)=%i \n',nNew,nOld,nSame,nNotSame,nNew-nSame)
    
    for iVar=1:nVar
        varargout{iVar}(same)=varargin{iVar}(ID(same));
    end
    
    
    
    %% Now check if any locations are remaining
    %
    
    if (nNew-nSame)>0
        
        % For the remaining (x2,y2) locations:
        % 1) if outside values have been defined, use those values
        % 2) elsewere use scattered interpolant
        %
        
        % Are any of the remaining locations within the old mesh?
        
        ID = pointLocation(MUA1.TR,[x2(~same) y2(~same)]) ;  % ID is NaN if outside all triangles, 
        NodesNotSame=find(~same);
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
        shift=1e-5;
        IDTest = pointLocation(MUA1.TR,[x2(NodesOutside)+shift y2(NodesOutside)+2*shift]) ;
        
        NewInsideAndNotSameNodes=NodesOutside(~isnan(IDTest));
        %  Add those nodes to the right set
        if ~isempty(NewInsideAndNotSameNodes)
            NodesInsideAndNotSame=[NodesInsideAndNotSame;NewInsideAndNotSameNodes];
            NodesOutside=setdiff(NodesOutside,NewInsideAndNotSameNodes);
        end
        
        IDTest = pointLocation(MUA1.TR,[x2(NodesOutside)-2*shift y2(NodesOutside)+shift]) ;
        NewInsideAndNotSameNodes=NodesOutside(~isnan(IDTest));
        %  Add those nodes to the right set
        if ~isempty(NewInsideAndNotSameNodes)
            NodesInsideAndNotSame=[NodesInsideAndNotSame;NewInsideAndNotSameNodes];
            NodesOutside=setdiff(NodesOutside,NewInsideAndNotSameNodes);
        end
        
        % And finally, are any of the outside nodes actually on the mesh boundary?
        NodesOnBoundary = DistanceToLineSegment([x2(NodesOutside) y2(NodesOutside)],[MUA1.Boundary.x MUA1.Boundary.y],[],1000*eps);
        
        %  Add any ouside nodes on bounday to the set of inside nodes
        if ~isempty(NodesOnBoundary)
            NodesInsideAndNotSame=[NodesInsideAndNotSame;NodesOutside(NodesOnBoundary)];
            NodesOutside(NodesOnBoundary)=[];
        end
        
        fprintf('#Outside=%i \t   #Inside and not same=%i \n',numel(NodesOutside),numel(NodesInsideAndNotSame))
        
        if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots
            FindOrCreateFigure("-Old and new nodes-");
            tt=axis;
            hold off
            PlotMuaMesh(CtrlVar,MUA1)
            if ~isequal(tt,[0 1 0 1])
                axis(tt)
            end
            hold on
            p2=plot(x2(NodesOutside)/CtrlVar.PlotXYscale,y2(NodesOutside)/CtrlVar.PlotXYscale,'ob');
            p3=plot(x2(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,y2(NodesInsideAndNotSame)/CtrlVar.PlotXYscale,'or');
            if ~isempty(p2) && ~isempty(p3)
                legend([p2 p3],'Outside','Inside')'northeastoutside'
            end
            axis tight
            hold off
        end
        
        F = scatteredInterpolant();
        F.Points=MUA1.coordinates;
        
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