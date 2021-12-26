function [xc,yc,LSF,xcEdges,ycEdges]=CalvingFrontLevelSetGeometricalInitialisation(CtrlVar,MUA,Xc,Yc,LSF,options)

% Takes initial calving front data (Xc,Yc) and a level set function LSF having correct sign, and finds the calving front
% locations (xcEdges,ycEdges) along the element edges of MUA.
%
% On input only the sign of LSF must be correct. This is, for nodes inside/upstream of calving fronts LSF must be positive,
% and for nodes outside/downstream of calving fronts, negative.
%
% On return LSF has be initialized geometrically based on the distance from nodes to the (xcEdges,ycEdges) intersections
% along element edges
%
% (xc,yc) are the calving front positions corresponding to the zero level of the level-set function
%
% Ideally, on  return (xc,yc) should be close to (Xc,Yc),  but figuresin out how to initialise the LSF in such a way that
% this is exactly correct has defeted me.
%
%
% The "EleEdges" method is slower than "InputPoints" method, and possibly in practice the resulting differences are not that
% significant.
%

arguments
    CtrlVar struct
    MUA     struct
    Xc      (:,1) double
    Yc      (:,1) double
    LSF     (:,1) double
    options.test logical = false
    options.plot logical = false
    options.method char {mustBeMember(options.method,{'InputPoints','EleEdges'})} = "EleEdges"
    options.CalvingFrontPointDistance (1,1) double =10e3; % this is the distance between the points defining the calving front
end


UserVar=[] ; RunInfo=[];

if nargin > 4 && ~isempty(Xc)
    P1=[Xc(:) Yc(:)] ;
end

if options.test

    load("TestingInterfaceCalculation.mat","MUA")
    P1=1000*[-500 -500 ; -500 500 ; 500 500  ; 500 -500 ;  ...
        200 -500 ;  200 -100 ; -200 -100 ; -200 -500 ;   ...
        -500 -500] ;
    CtrlVar.PlotXYscale=1000;

    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    inside = inpoly2([x y],P1);
    LSF=ones(MUA.Nnodes,1);
    LSF(~inside)=-1;

end

Nele=MUA.Nele;
con=MUA.connectivity;
coo=MUA.coordinates;

% To speed up the calculation:
%
% 1) Get rid of any points of the calving profile outside the mesh domain
io=inpoly2(P1, [MUA.Boundary.x MUA.Boundary.y]);
P1=P1(io,:) ;

% and 2) subsample the points along the calving profile
CtrlVar.GLtension=1 ; CtrlVar.GLds=options.CalvingFrontPointDistance;
[xsmooth,ysmooth]= Smooth2dPos(P1(:,1),P1(:,2),CtrlVar);
P1=[xsmooth(:) ysmooth(:)] ;

switch options.method

    case "EleEdges"

        %
        %  Here the intersections between individual ele edges and the calving front profile as defined by (Xc,Yc)
        %  are determined, and the LSF is then initialized as the distance to those intersection points

        EleEdges=[con con(:,1)]; % These are element edges, this is only for 3-node elements.
        % In general, these should be the nodal numbers of the corner nodes
        % so that this traces a closed loop for every element.

        P2=NaN(5*Nele,2);
        ii=1;
        for I=1:Nele

            P2(ii:ii+3,:)=coo(EleEdges(I,:),:);
            P2(ii+4,:)=[NaN NaN ];
            ii=ii+5 ;

        end
        P2(end,:)=[] ;

        % Now calculate the intersection points using this matlab function
        [xcEdges,ycEdges]=polyxpoly(P1(:,1),P1(:,2),P2(:,1),P2(:,2));

    case "InputPoints"

        % Here I simply take the input profile, so not much to do but getting rid
        % of points not needed

        xcEdges=P1(:,1); ycEdges=P1(:,2);

end


if options.plot
    LSFonInput=LSF;
    CtrlVar.PlotGLs=0;
    [xcOnInput,yconInput]=PlotCalvingFronts(CtrlVar,MUA,LSFonInput,color="k",LineStyle="--");
end

LSF=SignedDistUpdate(UserVar,RunInfo,CtrlVar,MUA,LSF,xcEdges,ycEdges);

CtrlVar.PlotGLs=0;  % to suppress the plot
[xc,yc]=PlotCalvingFronts(CtrlVar,MUA,LSF,color="k",LineStyle="--");

if options.plot

    FindOrCreateFigure("After: LSF and calving front")

    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,LSF/CtrlVar.PlotXYscale);
    hold on

    PlotMuaMesh(CtrlVar,MUA,[],"w");
    %tt=axis;
    plot(P1(:,1)/CtrlVar.PlotXYscale,P1(:,2)/CtrlVar.PlotXYscale,'-g.',LineWidth=2)
    hold on
    %axis(tt)
    plot(xc/CtrlVar.PlotXYscale,yc/CtrlVar.PlotXYscale,'-k.',LineWidth=2)
    axis equal

    FindOrCreateFigure("Before: LSF and calving front")

    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,LSFonInput/CtrlVar.PlotXYscale);
    hold on

    PlotMuaMesh(CtrlVar,MUA,[],"w");
    %tt=axis;
    plot(P1(:,1)/CtrlVar.PlotXYscale,P1(:,2)/CtrlVar.PlotXYscale,'-g.',LineWidth=2)
    hold on
    %axis(tt)
    plot(xcOnInput/CtrlVar.PlotXYscale,yconInput/CtrlVar.PlotXYscale,'-k.',LineWidth=2)
    axis equal



end

end