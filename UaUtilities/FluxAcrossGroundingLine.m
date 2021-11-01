function [qGL,qGLx,qGLy,Fub,Fvb,Fr,Fh,LakeNodes,GLgeo,ubGL,vbGL]=FluxAcrossGroundingLine(CtrlVar,MUA,GF,ub,vb,ud,vd,h,rho,Fub,Fvb,Fr,Fh,LakeNodes)

%  [qGL,qGLx,qGLy,Fub,Fvb,Fr,Fh,LakeNodes,GLgeo]=FluxAcrossGroundingLine(CtrlVar,MUA,GF,ub,vb,ud,vd,h,rho,Fub,Fvb,Fr,Fh,LakeNodes)
%
% Calculates flux across grounding line.
%
% Very simple routine. Use this as a template, you will most likely need to do some modification for your own purpose.
%
% Also, this routine is an utility routine and has not been tested and optimized in a similar way as typical Ua routines.
% 
% On return qGL is the vertically and horizontally integrated flux 
% over the depth h=s-b and the width dsGL=sqrt((GLgeo(:,3)-GLgeo(:,4)).^2+(GLgeo(:,5)-GLgeo(:,6)).^2);
% at xy-location (GLgeo(:,7),GLgeo(:,8)) 
%
% the units of flux are: velocity x thickness x width x density = m/s m m kg/m^3=kg/s (SI units)
%
%
% Fub, Fvb, Fr, Fh and LakeNodes are all optional inputs, and if not given are calculated and returned
% 
% If called several times for the same MUA, GF, ub, vb, ud, vd, h, and rho: 
% Give the Fub, Fvb, Fr, Fh, LakeNodes outputs from first call as inputs for all subsequent calls. 
% This has no effect on the results, but avoids recalculations and speeds things up
% 
%

if nargin<10
    Fub=[]; Fvb=[];Fr=[];Fh=[];LakeNodes=[];
end
    
CtrlVar.GLsubdivide=1;

if isempty(Fub)
    Fub=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ub);
    Fub.Method='natural' ; 
    Fvb=Fub; Fr=Fub; Fh=Fub;
    Fvb.Values=vb;
    Fr.Values=rho;
    Fh.Values=h;
%     Fvb=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),vb);
%     Fr=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),rho);
%     Fh=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),h);
    
end


% get rid of grounding lines around lakes
if nargin>13
    fprintf('Getting read of lake nodes (note: this is not a robust general approach, do check results.\n')
    if isempty(LakeNodes)
        [OceanNodes,LakeNodes]=LakeOrOcean(CtrlVar,GF,MUA.Boundary,MUA.connectivity,MUA.coordinates);
        GF.node(LakeNodes)=1;
    end
else
    LakeNodes=[] ; OceanNodes=[];
end

%GF.node(MUA.Boundary.Nodes)=0;


GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
xGL=GLgeo(:,7) ; yGL=GLgeo(:,8);
nxGL=GLgeo(:,9); nyGL=GLgeo(:,10);

rhoGL=Fr(xGL,yGL);
ubGL=Fub(xGL,yGL);
vbGL=Fvb(xGL,yGL);
hGL=Fh(xGL,yGL);
dsGL=sqrt((GLgeo(:,3)-GLgeo(:,4)).^2+(GLgeo(:,5)-GLgeo(:,6)).^2);
qGL=dsGL.*hGL.*(ubGL.*nxGL+vbGL.*nyGL).*rhoGL;  % units: m m (m/yr) * kg/m^3 =kg/yr
qGLx=qGL.*nxGL ; qGLy=qGL.*nyGL ;

% xa=GLgeo(:,3); xb=GLgeo(:,4); 
% ya=GLgeo(:,5); yb=GLgeo(:,6);
% qa=Fh(xa,ya).*Fr(xa,ya).*( Fub(xa,ya).*nxGL+Fvb(xa,ya).*nyGL); 
% qb=Fh(xb,yb).*Fr(xb,yb).*( Fub(xb,yb).*nxGL+Fvb(xb,yb).*nyGL); 
% qGL2=0.5*(qa+qb).*dsGL ; 
% qGL2x=qGL2.*nxGL ; qGL2y=qGL2.*nyGL ;

% 
% [ubGL2,vbGL2,rhoGL2,hGL2]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,MUA,xGL,yGL,ub,vb,rho,h);
% qGL2=dsGL.*hGL2.*(ubGL2.*nxGL+vbGL2.*nyGL).*rhoGL2;  % units: m m (m/yr) * kg/m^3 =kg/yr
% qGL2x=qGL2.*nxGL ; qGL2y=qGL2.*nyGL ;

return
%%


FindOrCreateFigure("GL Flux") ; 
PlotMuaMesh(CtrlVar,MUA) ; hold on ; 
plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'g','LineWidth',2); 
scale=1; hold on ;  quiver(GLgeo(:,7)/CtrlVar.PlotXYscale,GLgeo(:,8)/CtrlVar.PlotXYscale,qGLx,qGLy,scale,'color','r') ; 
axis equal

scatter(GLgeo(:,7)/CtrlVar.PlotXYscale,GLgeo(:,8)/CtrlVar.PlotXYscale,dsGL/CtrlVar.PlotXYscale) ; 



% 
% FindOrCreateFigure("GL Flux 2") ; 
% PlotMuaMesh(CtrlVar,MUA) ; hold on ; 
% plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'g','LineWidth',2); 
% scale=1; hold on ;  quiver(GLgeo(:,7)/CtrlVar.PlotXYscale,GLgeo(:,8)/CtrlVar.PlotXYscale,qGL2x,qGL2y,scale,'color','r') ; 
% axis equal
% 
% scatter(GLgeo(:,7)/CtrlVar.PlotXYscale,GLgeo(:,8)/CtrlVar.PlotXYscale,dsGL/CtrlVar.PlotXYscale) ; 
% %%
    
    
    
end