function [qGL,qGLx,qGLy,xEdge,yEdge,nxGL,nyGL,dsGL,xGLele,yGLele,iE,EN]=FluxAcrossEleGL(CtrlVar,MUA,GF,ub,vb,rho,h,DoPlots,iE,EN)

%%
% Calculates flux across grounding lines, with grounding lines defined as
% the free boundary (or free edges) of grounded elements.
%
%    [qGL,qGLx,qGLy,xEdge,yEdge,nxGL,nyGL,dsGL,xGLele,yGLele,iE,EN]=FluxAcrossEleGL(CtrlVar,MUA,GF,ub,vb,rho,h,DoPlots,iE,EN)
%
% Inputs:
% DoPlots : If true, some plots are created. (optional)
%  iE     : Corner nodes of grounding line edges. (optional)
%  EN     : Nodes of grounding line edges, including nodes along the edge.   (optional)
%
%  iE and EN can be left empty, in which case they will be calculated and
%  returned as outputs. To speed things up, use iE and EN from a previous call
%  if calling this routine again for the same MUA and GF.
%
% Outputs:
%  qGL            : vertically and horizontally integrated flux over each element edge defining the grounding line.
%  qGLx, qGLy     : the x, y components of a vector with the length qGL and normal to the grounding line (i.e. these are not x and y components of qGL)
%  xEdge, yEdge   : x and y coordinates of a center point of each GL edge.
%  nxGL, nyGL     : unit normal to GL edges
%  dsGL           : length of GL edges
%  xGLele, yGLele : ordered list of all corner points of GL edges, with NaN
%                   between individual GL segments. Useful for plotting GL.
%  iE             : Corner nodes of grounding line edges.
%  EN             : Nodes of grounding line edges, including nodes along the edge.
%
% Note: For 3-node elements, iE and EN are not needed.
%
% Note: qGL has the units: speed x ice-thickness x ice-density width = m/yr m kg/m^3 m = kg/yr
% (if using m, yr, kg as units for distance, time and mass.)
%
% Note: If using repeatedly for higher-order elements and the same MUA and GF,
%       give returned iE and EN from previous calls as inputs.
%
% Note: This routine uses a (slightly) different representation of the grounding line as that
% given by, for example, PlotGroundingLine.m which is based on interpolation of
% the GF nodal mask.
%
% Note: Currently only implemented for 3 and 6 node elements.

if nargin<8  || isempty(DoPlots)
    DoPlots=0;
end

if nargin<9
    iE=[] ;
    EN=[];
end


[xGLele,yGLele,triGR,FB,xEdge,yEdge,nxGL,nyGL,dsGL]=EleBasedGL(CtrlVar,MUA,GF,DoPlots);


switch MUA.nod
    
    case 3
        
        uhr=ub.*h.*rho;
        vhr=vb.*h.*rho;
        qxS=dsGL.*mean(uhr(FB),2) ;  % for a linear triangle element this is correct
        qyS=dsGL.*mean(vhr(FB),2) ;
        qGL=qxS.*nxGL+qyS.*nyGL; % units: m m (m/yr) * kg/m^3 =kg/yr
        qGLx=qGL.*nxGL ;
        qGLy=qGL.*nyGL ;
        EN=[FB(:,1) FB(:,2)];
    case 6
        
        % a bit more accurate estimate for 6-node elements. Here I used the two corner
        % nodes and the side node and integrate with Simpson's rule, which is second order
        % just like the form functions. So this should be pretty much exact.
        
        if nargin <10 || isempty(EN)
            if nargin<9 || isempty(iE)
                iE=FindEdgeNodesGivenTwoCornerNodes(MUA,FB);  % this can be slow
            end
            
            if size(iE,1)~=size(FB,1)
                fprintf('iE given as input does have a number of elements that is consistent with MUA, iE is therefore recalculated.\n')
                iE=FindEdgeNodesGivenTwoCornerNodes(MUA,FB);  % this can be slow
            end
            
            EN=[FB(:,1) iE(:) FB(:,2)];
            
        end
        
        uhr=ub.*h.*rho;
        vhr=vb.*h.*rho;
        qxS=dsGL.*(uhr(EN(:,1))+4*uhr(EN(:,2))+uhr(EN(:,3)))/6;  % Simpson's rule
        qyS=dsGL.*(vhr(EN(:,1))+4*vhr(EN(:,2))+vhr(EN(:,3)))/6;
        
        % I'm integrating the flux in x and y directions
        % and then projecting along normal and tangential directions.
        % This will only be correct for straight edges within each triangle
        %
        
        qGL=qxS.*nxGL+qyS.*nyGL;
        
        qGLx=qGL.*nxGL ;
        qGLy=qGL.*nyGL ;
        
    case 10
        
        error('Case for 10 node elements not yet implemented')
        
end


if DoPlots
    %%
    figure
    plot(xGLele,yGLele)
    hold on
    Par.VelPlotIntervalSpacing='log10';
    Par.RelativeVelArrowSize=10;
    Par.QuiverColorPowRange=2;
    % assuming units: m/yr, m, kg/m^3  for velocity, distance and density
    % q calculated has the units:  m/yr m m kg/m^3 = kg/yr
    % To get Gt divide by 1e9
    cbar=QuiverColorGHG(xEdge,yEdge,qGLx/1e9,qGLy/1e9,Par);
    title(cbar,'kg/yr')   ;
    axis equal
    
end

end
