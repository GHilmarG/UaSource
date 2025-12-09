function [Tri]=PlotNodalVariableAsTriSurface(CtrlVar,MUA,Tri,Variable,AspectRatio)

% [Tri]=PlotNodalVariableAsTriSurface(CtrlVar,MUA,Tri,Variable,AspectRatio)
%
% AspectRatio is optional
% Tri can be left empty
%
%  Example
%
% figure ; PlotNodalVariableAsTriSurface(CtrlVar,MUA,[],log10(C),1000);
%
%

if isempty(Tri)
    [TRI,DT]=CreateTRI(MUA);
end

if isempty(CtrlVar)
    CtrlVar.PlotXYscale=1;
end

if nargin<5 || isempty(AspectRatio)
    % if no aspect ratio given, try to come up with something sensible
    if max(Variable) > min(Variable)
        AspectRatio=1e-2*(max(MUA.coordinates(:,1))-min(MUA.coordinates(:,1)))/(max(Variable)-min(Variable));
    else
        AspectRatio=1;
    end
end

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;

trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,Variable,'EdgeColor','none') ;
axis equal ; tt=daspect ;
daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)*CtrlVar.PlotXYscale/10/AspectRatio]);
axis tight
lighting phong ;

end