function lgd=PlotBoundaryConditions(CtrlVar,MUA,BCs,varargin)
%%
%  PlotBoundaryConditions(CtrlVar,MUA,BCs,varargin)
%
%  Gives a graphical representation of boundary conditions. 
%
%  varargin is passed on to PlotMuaMesh
%
%   lgd :    legend object
%
% Example:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   figure ; PlotBoundaryConditions(CtrlVar,MUA,BCs,'k')
%
%%



CtrlVar.PlotMesh=1; CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
I=0; L=[];
PlotMuaMesh(CtrlVar,MUA,[],varargin{:})  ; hold on

x=MUA.coordinates(:,1)/CtrlVar.PlotXYscale; y=MUA.coordinates(:,2)/CtrlVar.PlotXYscale;

velscale=min([max(x)-min(x) ; max(y)-min(y)])/30;
velscale=velscale*CtrlVar.BoundaryConditionsFixedNodeArrowScale;
headscale=0.3; sharp=0.3; head=1; lw=1; io=-1; col='r';

if strcmpi(CtrlVar.FlowApproximation,'SSTREAM') || strcmpi(CtrlVar.FlowApproximation,'Hybrid')
    
    [nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
    
    if ~isempty(BCs.ubFixedNode)
        
        % I want the arrows to point towards the nodes and be outside of the model
        % domain. Hence I check what the nodal normal is.
        xfixed=x(BCs.ubFixedNode); yfixed=y(BCs.ubFixedNode);
        xNorm=Nx(BCs.ubFixedNode)./abs(Nx(BCs.ubFixedNode)); 
        xNorm(isnan(xNorm))=1; % if the normal is exactly at 90 degrees to the x-direction then Nx will be zero and I get NaN. In this case just set it to either +1 or -1.
        yNorm=xfixed*0;
        ghg_arrow(xfixed,yfixed,xNorm,yNorm,velscale,headscale,sharp,head,col,lw,io);

    end
   
    if ~isempty(BCs.vbFixedNode)

        xfixed=x(BCs.vbFixedNode); yfixed=y(BCs.vbFixedNode);
        xNorm=xfixed*0;
        yNorm=Ny(BCs.vbFixedNode)./abs(Ny(BCs.vbFixedNode)); 
        yNorm(isnan(yNorm))=1;
        ghg_arrow(xfixed,yfixed,xNorm,yNorm,velscale,headscale,sharp,head,col,lw,io);
        %ghg_arrow(xfixed,yfixed,xfixed*0,yfixed*0+1,velscale,headscale,sharp,head,col,lw,io);

    end
    
    if ~isempty(BCs.ubvbFixedNormalNode)
        xfixed=x(BCs.ubvbFixedNormalNode); yfixed=y(BCs.ubvbFixedNormalNode);
        
        col='c';
        ghg_arrow(xfixed,yfixed,Nx(BCs.ubvbFixedNormalNode),Ny(BCs.ubvbFixedNormalNode),velscale,headscale,sharp,head,col,lw,io);
    end
    % plot ties
    
    if ~isempty(BCs.ubTiedNodeA)
        for I=1:numel(BCs.ubTiedNodeA)
            plot(x(BCs.ubTiedNodeA(I)),y(BCs.ubTiedNodeA(I)),'ob')
            plot(x(BCs.ubTiedNodeB(I)),y(BCs.ubTiedNodeB(I)),'xb')
            plot([x(BCs.ubTiedNodeA(I))  x(BCs.ubTiedNodeB(I))],[y(BCs.ubTiedNodeA(I))  y(BCs.ubTiedNodeB(I))],'b--')
        end
    end
    
    if ~isempty(BCs.vbTiedNodeA)
        for I=1:numel(BCs.vbTiedNodeA)
            plot(x(BCs.vbTiedNodeA(I)),y(BCs.vbTiedNodeA(I)),'xr')
            plot(x(BCs.vbTiedNodeB(I)),y(BCs.vbTiedNodeB(I)),'^r')
            plot([x(BCs.vbTiedNodeA(I))  x(BCs.vbTiedNodeB(I))],[y(BCs.vbTiedNodeA(I))  y(BCs.vbTiedNodeB(I))],'r-.')
        end
    end
    
end

if strcmp(CtrlVar.FlowApproximation,'SSHEET')
    
    if ~isempty(BCs.udFixedNode)
        xfixed=x(BCs.udFixedNode); yfixed=y(BCs.udFixedNode);
        ghg_arrow(xfixed,yfixed,xfixed*0+1,yfixed*0,velscale,headscale,sharp,head,col,lw,io);
    end
    
    if ~isempty(BCs.vdFixedNode)
        xfixed=x(BCs.vdFixedNode); yfixed=y(BCs.vdFixedNode);
        ghg_arrow(xfixed,yfixed,xfixed*0,yfixed*0+1,velscale,headscale,sharp,head,col,lw,io);
    end
    
    if ~isempty(BCs.udvdFixedNormalNode)
        xfixed=x(BCs.udvdFixedNormalNode); yfixed=y(BCs.udvdFixedNormalNode);
        [nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
        col='c';
        ghg_arrow(xfixed,yfixed,Nx(BCs.udvdFixedNormalNode),Ny(BCs.udvdFixedNormalNode),velscale,headscale,sharp,head,col,lw,io);
    end
    % plot ties
    
    if ~isempty(BCs.ubTiedNodeA)
        for I=1:numel(BCs.vbTiedNodeA)
            plot(x(BCs.ubTiedNodeA(I)),y(BCs.ubTiedNodeA(I)),'ob')
            plot(x(BCs.ubTiedNodeB(I)),y(BCs.ubTiedNodeB(I)),'xb')
            plot([x(BCs.ubTiedNodeA(I))  x(BCs.ubTiedNodeB(I))],[y(BCs.ubTiedNodeA(I))  y(BCs.ubTiedNodeB(I))],'b--')
        end
    end
    
    if ~isempty(BCs.vbTiedNodeA)
        for I=1:numel(BCs.vbTiedNodeA)
            plot(x(BCs.vbTiedNodeA(I)),y(BCs.vbTiedNodeA(I)),'xr')
            plot(x(BCs.vbTiedNodeB(I)),y(BCs.vbTiedNodeB(I)),'^r')
            plot([x(BCs.vbTiedNodeA(I))  x(BCs.vbTiedNodeB(I))],[y(BCs.vbTiedNodeA(I))  y(BCs.vbTiedNodeB(I))],'r-.')
        end
    end
    
end

if ~isempty(BCs.hTiedNodeA)
    for I=1:numel(BCs.vbTiedNodeA)
        plot(x(BCs.hTiedNodeA(I)),y(BCs.hTiedNodeA(I)),'sg')
        plot(x(BCs.hTiedNodeB(I)),y(BCs.hTiedNodeB(I)),'dg')
        plot([x(BCs.hTiedNodeA(I))  x(BCs.hTiedNodeB(I))], [y(BCs.hTiedNodeA(I))  y(BCs.hTiedNodeB(I))],'g:')
    end
end

if ~isempty(BCs.hFixedNode)
    xfixed=x(BCs.hFixedNode); yfixed=y(BCs.hFixedNode);
    plot(xfixed,yfixed,'oc','MarkerFaceColor','c')
end

if ~isempty(BCs.hPosNode)
    xfixed=x(BCs.hPosNode); yfixed=y(BCs.hPosNode);
    I=I+1 ; L{I}=plot(xfixed,yfixed,'*b','MarkerFaceColor','b','DisplayName','$h$'); 
end


if ~isempty(BCs.LSFFixedNode)
    xfixed=x(BCs.LSFFixedNode); yfixed=y(BCs.LSFFixedNode);
    I=I+1; L{I}=plot(xfixed,yfixed,'or','MarkerFaceColor','b','DisplayName','$\varphi$'); 
end




if strcmpi(CtrlVar.FlowApproximation,'SSTREAM') || strcmpi(CtrlVar.FlowApproximation,'Hybrid')
title(...
    sprintf('Boundary conditions: \n Arrows represent fixed ub,vb, and normal velocites (%i,%i,%i). \n Cyan and blue symbols show where the thickness is prescribed/constrained (%i,%i) \n Blue, red and grean lines are (ub,vb,h) nodal ties (%i,%i,%i)',...
    numel(BCs.ubFixedNode),numel(BCs.vbFixedNode),numel(BCs.ubvbFixedNormalNode),numel(BCs.hFixedNode),numel(BCs.hPosNode),numel(BCs.ubTiedNodeA),numel(BCs.vbTiedNodeA),numel(BCs.hTiedNodeA)),...
    'FontSize',9)
elseif strcmp(CtrlVar.FlowApproximation,'SSHEET')
    title(...
    sprintf('Boundary conditions: \n Arrows represent fixed ud,vd, and normal velocites (%i,%i,%i). \n Cyan and blue symbols show where the thickness is prescribed/constrained (%i,%i) \n Blue, red and grean lines are (ud,vd,h) nodal ties (%i,%i,%i)',...
    numel(BCs.udFixedNode),numel(BCs.vdFixedNode),numel(BCs.udvdFixedNormalNode),numel(BCs.hFixedNode),numel(BCs.hPosNode),numel(BCs.ubTiedNodeA),numel(BCs.vbTiedNodeA),numel(BCs.hTiedNodeA)),...
    'FontSize',9)
end

if ~isempty(L)
    lgd=legend([L{:}],'interpreter','latex');
else
    lgd=[]; 
end


end
