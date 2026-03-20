


function [B,UserVar,CtrlVar,MUA,F]=AnalyticalBedgeometries(UserVar,CtrlVar,MUA,F,options)

% F can be a structure or [x y] array
% function [B,UserVar,CtrlVar,MUA,F]=AnalyticalBedgeometries(UserVar,CtrlVar,MUA,[x y],options)


arguments
    UserVar struct=[]
    CtrlVar struct=[]
    MUA     struct=[]
    F       {mustBeA(F,{'struct','UaFields','numeric'})}=[]
    options.Plot  logical = false
    options.BedName string="-Thule-"
end



if nargin == 0 || isempty(MUA)  || isempty(F)

    %%
    options.Plot=true;
    R=1000e3 ;
    theta=linspace(0,2*pi,100);
    x=R*cos(theta); y=R*sin(theta) ; x(end)=[] ; y(end )=[];


    UserVar=[];
    CtrlVar=Ua2D_DefaultParameters(); %
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=1;
    % Note; When creating this mesh using Úa, only the following
    % three lines are required in the Ua2D_InitialUserInput.m
    CtrlVar.MeshSizeMax=50e3;
    CtrlVar.MeshSizeMin=10e3;
    CtrlVar.MeshSize=10e3 ;

    MeshBoundaryCoordinates=[x(:) y(:)];

    CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
    % Now generate mesh (When using Úa this is done internally, no such call
    % then needed).


    [UserVar,MUA]=genmesh2d(UserVar,CtrlVar);
    figure ; PlotMuaMesh(CtrlVar,MUA); drawnow


    F.x= MUA.coordinates(:,1); F.y= MUA.coordinates(:,2);



end

if isa(F,"struct")

    x=F.x ; y=F.y ;
  
else

    x=F(:,1);
    y=F(:,2);

end

B=Bfunc(x,y,options.BedName) ;

if options.Plot


    FindOrCreateFigure("B Thule") ;
    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,B);
    xlabel("x (km)",Interpreter="latex") ;
    ylabel("y (km)",Interpreter="latex")
    colormap(othercolor("Mdarkterrain",25))
    axis tight
    title(cbar,"(m a.s.l.)",Interpreter="latex")

    %ax=gca; exportgraphics(ax,'ThuleB.pdf')

    x=linspace(0,1000e3,100) ; y=x*0;
    Bxprofile=Bfunc(x,y,options.BedName);

    FindOrCreateFigure("B: x-profile") ;
    plot(x/1000,Bxprofile) ;
    xlabel("x (km)") ; ylabel("B (m)")

    y=linspace(0,1000e3,100) ; x=y*0;
    Bxprofile=Bfunc(x,y,options.BedName);

    FindOrCreateFigure("B: y-profile") ;
    plot(y/1000,Bxprofile) ;
    xlabel("y (km)") ; ylabel("B (m)")


end

if ~nargout   % A trick to suppress any function output if no output requested. No need to suppress output using ;
    clearvars B
end

end



%% local function

function B=Bfunc(x,y,BedName)





%%   new parameter set
R=800e3 ;
Bc=900;
Bl=-2000;
Ba=1100;
%%

switch BedName

    case "-Thule-"

        r=sqrt(x.*x+y.*y) ;
        theta=atan2(y,x);


        rc=0;
        l=R -  cos(2*theta).*R/2 ;            % theta-dependent wavelength
        a=Bc - (Bc-Bl)*(r-rc).^2./(R-rc).^2;  % quadratic term in r
        B=Ba*cos(3*pi*r./l)+a ;               %

    case "-ThuleNS-"

        r=sqrt(x.*x+y.*y) ;
        theta=pi/2 ;

        rc=0;
        l=R -  cos(2*theta).*R/2 ;            % theta-dependent wavelength
        a=Bc - (Bc-Bl)*(r-rc).^2./(R-rc).^2;  % quadratic term in r
        B=Ba*cos(3*pi*r./l)+a ;

    case "-MismBed-"

        B=MismBed(x,y);

    otherwise

        error("Case not found")

end

end




