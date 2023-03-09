function [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,varargin)



%% Plots grounding lines
%
% To plot grounding lines over FE mesh based on the floating mask GL:
%
%  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,varargin)
%
%
% To plot (most of the) grounding lines based on the Bedmachine data:
%
%  [cGL,yGL]=PlotGroundingLines(CtrlVarInRestartFile,"Bedmachine");
%
%  [cGL,yGL]=PlotGroundingLines([],"Bedmachine");    
%
% When plotting grounding lines over the mesh defined by MUA and based on GF, the only required inputs are:
%
%   MUA GF
%
% Other fields can be left empty. However, if the grounding line needs to be plotted repeatedly for same MUA and GF, entering
% GLgeo, xGL and yGL, obtained as outputs from a previouis call, will speed things up.
%
% varargin is passed over the the plot function and can be any input
% accepted by the matlap plot function.
%
% *Examples:*
%
% Plot grounding lines with minimum of required input:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','GF','CtrlVar')
%   Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity); Tlength=sqrt(2*Tarea) ;
%   figure ;
%   [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF) ;
%
% Plot grounding lines in red over the computational mesh in black:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','GF','CtrlVar')
%   figure
%   CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
%   CtrlVar.PlotIndividualGLs=1 ;
%   PlotMuaMesh(CtrlVar,MUA,[],'k') ;
%   hold on
%   PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'color','r','LineWidth',2);
%
% Plot grounding lines twice using outputs from first call in the second call:
%
%    load('MUA-PIG-TWG-Example.mat','MUA','BCs','GF','CtrlVar')
%    GLgeo=[] ; xGL=[], yGL=[];
%    figure
%    tic; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL); toc
%    tic; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL); toc
%
% To plot grounding lines for higher-order elements with sub-element resolution
% set CtrlVar.GLsubdivide=1. Doing so causes higher-order elements to be split
% into sub-elements allowing for more accurate representation of the grounding
% line within an element.
%
% If xGL and yGL are empty on input, they are calculated from MUA and GLgeo.
%
% If GLgeo, xGL and yGL are all empty on input, then GLgeo is first calculated
% from from GF and MUA, and then xGL and yGL are calculated from GLgeo.
%
% If GLgeo, xGL and yGL are all empty, then MUA and GF can not be empty.
%
% If xGL and yGL are not empty, then all other fields can be left empty.
%
%
% Once xGL and yGL have been calculated, one can also simply plot the grounding lines as:
%
%   plot(xGL,yGL)
%
% To plot individual grounding lines in different colours set
%
%   CtrlVar.PlotIndividualGLs=1
%
% To just calculate the grounding line, but not plot it set
%
%   CtrlVar.PlotGLs=0;
%   [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
%
% Plot grounding lines based on the Bedmachine mask
%
%   
%
%
% See also: PlotMuaBoundary, EleBasedGL
%

narginchk(0,inf)


if nargin==0 

    CtrlVar=[]; 
    MUA="Bedmachine" ; 

end

if isempty(CtrlVar)
    CtrlVar.PlotXYscale=1000;
    CtrlVar.PlotIndividualGLs=0;
    CtrlVar.PlotGLs=1;
end

if ~isfield(CtrlVar,"PlotGLs")
    CtrlVar.PlotGLs=1;
end

if ~isfield(CtrlVar,"PlotIndividualGLs")
    CtrlVar.PlotIndividualGLs=0;
end



if isstring(MUA)

    GroundingLineDataSet=MUA ;

    switch GroundingLineDataSet

        case "Bedmachine"

            load("GroundingLineForAntarcticaBasedOnBedmachine.mat","xGL","yGL") ;

        case "Bindschadler"

            load('GroundingLinesOfAntarticaFromBobBindschadler','xGL','yGL')  ;

        otherwise

            error("case not found")

    end

    if CtrlVar.PlotGLs

        tt=axis;
        plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,varargin{:}) ;
        ax=gca; ax.DataAspectRatio=[1 1 1];

        if ~isequal(tt,[0 1 0 1])
            axis(tt)
        end
    end

    GLgeo=[] ;
    return

end



if isempty(GF)

    return
end

if ~isfield(CtrlVar,'PlotXYscale') ; CtrlVar.PlotXYscale=1 ; end
if ~isfield(CtrlVar,'PlotGLs') ; CtrlVar.PlotGLs=1 ; end
if ~isfield(CtrlVar,'PlotIndividualGLs') ; CtrlVar.PlotIndividualGLs=0 ; end
if ~isfield(CtrlVar,'LineUpGLs') ; CtrlVar.LineUpGLs=1; end


 if CtrlVar.PlotGLs  % if plotting, always line up grounding lines
     CtrlVar.LineUpGLs=1;
 end

if nargin<4 || isempty(GLgeo)

    GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);

end

if nargin<6 || ( isempty(xGL) || isempty(yGL))


    if CtrlVar.LineUpGLs
        xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ; ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
        [xGL,yGL]=LineUpEdges2(CtrlVar,xa,xb,ya,yb);

        %% get rid of duplicats and almost duplicates

    else
        xGL=[GLgeo(:,3)  ; GLgeo(:,4) ] ;
        yGL=[GLgeo(:,5)  ; GLgeo(:,6) ] ;
       %  temp=unique([xGL yGL],'rows') ;
        temp=uniquetol([xGL yGL],1000*eps,ByRows=true) ;

        xGL=temp(:,1) ;  yGL=temp(:,2) ;
    end

end



if CtrlVar.PlotGLs

    if ~CtrlVar.PlotIndividualGLs

        plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,varargin{:}) ;
        ax=gca; ax.DataAspectRatio=[1 1 1];

    else

        i=0;
        I=find(isnan(xGL)) ;
        I=[1;I(:)];
        col=['b','r','c','g','k','m'];
        %col=['b','r'];
        for ii=1:numel(I)-1
            i=i+1;

            plot(xGL(I(ii):I(ii+1))/CtrlVar.PlotXYscale,yGL(I(ii):I(ii+1))/CtrlVar.PlotXYscale,'color',col(i),varargin{:}) ;
            axis equal ; hold on ;
            if i==numel(col) ; i=0 ; end
        end
        ax=gca; ax.DataAspectRatio=[1 1 1];
    end
end



end