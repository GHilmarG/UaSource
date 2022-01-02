function [xc,yc,LSF,xcEdges,ycEdges,ShapeDifference]=CalvingFrontLevelSetGeometricalInitialisation(CtrlVar,MUA,Xc,Yc,LSF,options)

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
% Ideally, on  return (xc,yc) should be close to (Xc,Yc)
%
% The "EleEdges" method is slower than "InputPoints" method, and possibly in practice the resulting differences are not that
% significant.
%
% When using the "InputPoints" method, it might be good to resample the calving front at equal intervals by setting
%
%
%   ResampleCalvingFront=true
%   CalvingFrontPointDistance=1000 ; % the distance between points along the re-sampled calving front
%
%


arguments
    CtrlVar struct
    MUA     struct
    Xc      (:,1) double
    Yc      (:,1) double
    LSF     (:,1) double
    options.test logical = false
    options.plot logical = false
    options.method char {mustBeMember(options.method,{'InputPoints','EleEdges'})} = "InputPoints"
    options.ResampleCalvingFront logical = false
    options.CalvingFrontPointDistance (1,1) double = 1e3; % this is the (default) distance between the points defining the calving front
    options.GetRidOfCalvingFrontOutsideComputationalDomain logical = false
end

ShapeDifference=[];


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

% To speed up the calculation consider to:
% 1) get rid of any sections of the calving front outside of the computational boundary
% 2) subsample points along the calving front.
%
% But presumably all such manipulations with the input data are better done ahead of
% the call.

if options.GetRidOfCalvingFrontOutsideComputationalDomain
    % 1) Get rid of any points of the calving profile outside the mesh domain
    % But be carfulle with this as it does not deal with more than one calving front!
    io=inpoly2(P1, [MUA.Boundary.x MUA.Boundary.y]);
    P1=P1(io,:) ;
end

if options.ResampleCalvingFront
    % and 2) subsample the points along the calving profile

    ds=options.CalvingFrontPointDistance;

    nPoints=size(P1,1);

    % Check if there are several seperate calving fronts, and resample each individually
    iNaN=find(isnan(P1(:,1)));
    ii=[0;iNaN;nPoints+1] ;

    xi=[]; yi=[];
    for k=1:numel(ii)-1

        px=P1(ii(k)+1:ii(k+1)-1,1);
        py=P1(ii(k)+1:ii(k+1)-1,2);


        arclength = sum(sqrt(sum(diff([px py],[],1).^2,2)));
        Npoints=round(arclength/ds);

        pt= interparc(Npoints,px,py,'linear') ;
        xi=[xi;nan;pt(:,1)];
        yi=[yi;nan;pt(:,2)];
    end

    xi(end)=[]; yi(end)=[];

    P1=[xi(:) yi(:)] ;


end




switch options.method

    case "EleEdges"

        %
        %  Here the intersections between individual ele edges and the calving front profile as defined by (Xc,Yc)
        %  are determined, and the LSF is then initialized as the distance to those intersection points

        switch MUA.nod

            case 3

                EleEdges=con(:,[1 2 3 1]); % These are element edges, this is only for 3-node elements.
                % In general, these should be the nodal numbers of the corner nodes
                % so that this traces a closed loop for every element.

            case 6

                EleEdges=con(:,[1 3 5 1]); % These are element edges

                % This does not make it any more accurate, only slower
                %con=TriFE(MUA.connectivity);
                %EleEdges=con(:,[1 2 3 1]); % These are element edges, this is only for 3-node elements.



            case 10

                EleEdges=con(:,[1 4 7 1]); % These are element edges

                % This does not make it any more accurate, only slower
                %con=TriFE(MUA.connectivity);
                %EleEdges=con(:,[1 2 3 1]); % These are element edges, this is only for 3-node elements.

            otherwise

                error('not implemented')
        end

        nTri=size(con,1);
        P2=NaN(5*nTri,2);
        ii=1;
        for I=1:nTri

            P2(ii:ii+3,:)=coo(EleEdges(I,:),:);
            P2(ii+4,:)=[NaN NaN ];
            ii=ii+5 ;

        end
        P2(end,:)=[] ;

        % Now calculate the intersection points using this matlab function
        [xcEdges,ycEdges]=polyxpoly(P1(:,1),P1(:,2),P2(:,1),P2(:,2));

    case "InputPoints"

        % Here I simply take the input calving profile.
        %

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
    plot(Xc/CtrlVar.PlotXYscale,Yc/CtrlVar.PlotXYscale,'-g.',LineWidth=1,MarkerSize=12)
    hold on
    %axis(tt)
    plot(xc/CtrlVar.PlotXYscale,yc/CtrlVar.PlotXYscale,'-k.',LineWidth=1,MarkerSize=12)
    plot(xcEdges/CtrlVar.PlotXYscale,ycEdges/CtrlVar.PlotXYscale,'r.',LineWidth=1,MarkerSize=12)
    axis equal

    FindOrCreateFigure("Before: LSF and calving front")

    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,LSFonInput/CtrlVar.PlotXYscale);
    hold on

    PlotMuaMesh(CtrlVar,MUA,[],"w");
    %tt=axis;
    plot(Xc/CtrlVar.PlotXYscale,Yc/CtrlVar.PlotXYscale,'-g.',LineWidth=1)
    hold on
    %axis(tt)
    plot(xcOnInput/CtrlVar.PlotXYscale,yconInput/CtrlVar.PlotXYscale,'-k.',LineWidth=1,MarkerSize=12)


    axis equal

    if nargout>=6

        % Evaluate distance between prescribed calving front [Xc Yc] and the resulting calving front [xc yc], geometries.
        xy=[Xc Yc ; nan nan ; xc yc];
        ShapeDifference=polyshape(xy);
        FindOrCreateFigure("Shape Difference 2") ;
        plot(ShapeDifference) ; axis equal

    end

end

end