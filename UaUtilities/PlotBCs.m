
function PlotBCs(CtrlVar,coordinates,connectivity,Boundary,...
        ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,hfixednode,hfixedvalue,htiedA,htiedB,FixedNormalVelocityNode,FixedNormalVelocityValue)
    
 
   error('do not use. use PlotBoundaryCondtions instead')
    
    
    CtrlVar.PlotMesh=1;
 
    PlotFEmesh(coordinates,connectivity,CtrlVar); hold on
	
    coordinates=coordinates/CtrlVar.PlotXYscale;
    x=coordinates(:,1); y=coordinates(:,2);
    
	velscale=max([max(coordinates(:,1))-min(coordinates(:,1)) ; max(coordinates(:,2))-min(coordinates(:,2))])/20;
    velscale=velscale*CtrlVar.BoundaryConditionsFixedNodeArrowScale;
	headscale=0.3; sharp=0.3; head=1; lw=1; io=-1; col='r';
	
    if ~isempty(ufixednode)
        xfixed=coordinates(ufixednode,1); yfixed=coordinates(ufixednode,2);
        ghg_arrow(xfixed,yfixed,xfixed*0+1,yfixed*0,velscale,headscale,sharp,head,col,lw,io);
    end
    
    if ~isempty(vfixednode)
        xfixed=coordinates(vfixednode,1); yfixed=coordinates(vfixednode,2);
        ghg_arrow(xfixed,yfixed,xfixed*0,yfixed*0+1,velscale,headscale,sharp,head,col,lw,io);
    end
    
    if ~isempty(FixedNormalVelocityNode)
        xfixed=coordinates(FixedNormalVelocityNode,1); yfixed=coordinates(FixedNormalVelocityNode,2);
        [nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(connectivity,coordinates,Boundary.Edges);
        col='c';
        ghg_arrow(xfixed,yfixed,Nx(FixedNormalVelocityNode),Ny(FixedNormalVelocityNode),velscale,headscale,sharp,head,col,lw,io);
    end
    % plot ties
    
    if ~isempty(utiedA)
        for I=1:numel(vtiedA)
            plot(x(utiedA(I)),y(utiedA(I)),'ob')
            plot(x(utiedB(I)),y(utiedB(I)),'xb')
            plot([x(utiedA(I))  x(utiedB(I))],[y(utiedA(I))  y(utiedB(I))],'b--')
        end
    end
    
    if ~isempty(vtiedA)
        for I=1:numel(vtiedA)
            plot(x(vtiedA(I)),y(vtiedA(I)),'xr')  
            plot(x(vtiedB(I)),y(vtiedB(I)),'^r')  
            plot([x(vtiedA(I))  x(vtiedB(I))],[y(vtiedA(I))  y(vtiedB(I))],'r-.')
        end
    end
    
     if ~isempty(htiedA)
        for I=1:numel(vtiedA)
            plot(x(htiedA(I)),y(htiedA(I)),'sg')  
            plot(x(htiedB(I)),y(htiedB(I)),'dg')  
            plot([x(htiedA(I))  x(htiedB(I))], [y(htiedA(I))  y(htiedB(I))],'g:')
        end
    end
    
    if isfield(Boundary,'ElementsBCu')
        
        
        fe=connectivity(Boundary.ElementsBCu{1},Boundary.Edge{1})';
        plot(x(fe)+2*velscale/1e3, y(fe)+2*velscale/1e3, '-g', 'LineWidth',2) ;
        
        fe=connectivity(Boundary.ElementsBCu{2},Boundary.Edge{2})';
        plot(x(fe)+2*velscale/1e3, y(fe)+2*velscale/1e3, '-r', 'LineWidth',2) ;
        
        fe=connectivity(Boundary.ElementsBCu{3},Boundary.Edge{3})';
        plot(x(fe)+2*velscale/1e3, y(fe)+2*velscale/1e3, '-c', 'LineWidth',2) ;
        
        fe=connectivity(Boundary.ElementsBCv{1},Boundary.Edge{1})';
        plot(x(fe)+2*velscale/1e3, y(fe)-2*velscale/1e3, '-.g', 'LineWidth',2) ;
        
        fe=connectivity(Boundary.ElementsBCv{2},Boundary.Edge{2})';
        plot(x(fe)+2*velscale/1e3, y(fe)-2*velscale/1e3, '-.r', 'LineWidth',2) ;
        
        fe=connectivity(Boundary.ElementsBCv{3},Boundary.Edge{3})';
        plot(x(fe)+2*velscale/1e3, y(fe)-2*velscale/1e3, '-.c', 'LineWidth',2)
        
    end
    
    % plot edges
    %figure
    %fe=connectivity(Boundary.Elements{1},Boundary.Edge{1})'; plot(x(fe)+2*velscale/1e3, y(fe)+2*velscale/1e3, '-g', 'LineWidth',2) ; hold on
    %fe=connectivity(Boundary.Elements{2},Boundary.Edge{2})'; plot(x(fe)+2*velscale/1e3, y(fe)+2*velscale/1e3, '-r', 'LineWidth',2) ;
    %fe=connectivity(Boundary.Elements{3},Boundary.Edge{3})'; plot(x(fe)+2*velscale/1e3, y(fe)+2*velscale/1e3, '-c', 'LineWidth',2) ;
    %%


   if ~isempty(hfixednode)
        xfixed=coordinates(hfixednode,1); yfixed=coordinates(hfixednode,2);
        plot(xfixed,yfixed,'oc','MarkerFaceColor','c')
   end
   
   
   title(...
       sprintf('Boundary conditions: \n Arrows represent fixed u,v, and normal velocites (%i,%i,%i). \n Cyan filled circles show where the thickness is fixed (%i) \n Blue, red, and green lines are (u,v,h) nodal ties (%i,%i,%i)',...
       numel(ufixednode),numel(vfixednode),numel(FixedNormalVelocityNode),numel(hfixednode),numel(utiedA),numel(vtiedA),numel(htiedA)),...
       'FontSize',9)


end
