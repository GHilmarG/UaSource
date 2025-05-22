






function [cbar,xGL,yGL,xCF,yCF,CtrlVar,lg]=UaPlots(CtrlVar,MUA,F,Variable,options)

%%
%
% cbar=UaPlots(CtrlVar,MUA,F,Variable,options)
%
% Simple plot utility to plot variables and calving fronts and grounding lines as well.
%
% Note: Sometimes the default labels on the plots assume some typical
%       physical dimensions such as m/yr for velocities, and kPa for stresses.
% 
%
% Note:  To produce two velocity plots with the same scaling, use the CtrlVar from previous call again, but in the second
% call set
% 
%   CtrlVar.QuiverSameVelocityScalingsAsBefore=true;
%
%
% Returns grounding lines (xGL,yGL) and calving fronts (xCF,yCF).
%
%
%
% Examples:
%
%   UaPlots(CtrlVar,MUA,F,F.h)
%
%   UaPlots(CtrlVar,MUA,F,"-speed-")
%
%   UaPlots(CtrlVar,MUA,F,"-ubvb-")
%
%   cbar=UaPlots(CtrlVar,MUA,F,F.h,CalvingFrontColor="b",GroundingLineColor="k",GetRidOfValuesDownStreamOfCalvingFronts=false,ColorMap=jet);
%   title(cbar,"h (m)")
%
%   UaPlots(CtrlVar,MUA,F,"-log10speed-",CalvingFrontColor="b",GroundingLineColor="k",GetRidOfValuesDownStreamOfCalvingFronts=false,ColorMap=othercolor("YlGnBu8",100),PlotUnderMesh=true) ;
%
%
%   figetaInt=FindOrCreateFigure("eta Int") ; clf(figetaInt) ;
%   UaPlots(CtrlVar,MUA,F,"eta int",GetRidOfValuesDownStreamOfCalvingFronts=false) ;
%
%
% Plotting velocities other than those in F:
%
%     dub=F1.ub-F0.ub ; dvb=F1.vb-F0.vb ; 
%     UaPlots(CtrlVar,MUA,F1,[dub dvb],FigureTitle="(duv)")
%
% Log color scale:
% 
%    cbar=UaPlots(CtrlVar,MUA,F,abs(F.ab)); set(gca,'ColorScale','log') 
%
% Basal melt distribution, using log scale:
%
%    UaPlots(CtrlVar,MUA,F,"-log(ab)-")
%
%
% To plot velocities and set the range
%
%   CtrlVar.QuiverColorSpeedLimits=[0 2000];
%   UaPlots(CtrlVar,MUA,F,"-uv-",FigureTitle="velocities")
%
% Plot locations of min ice thickness:
%
%   UaPlots(CtrlVar,MUA,F,F.h,GetRidOfValuesDownStreamOfCalvingFronts=false,logColorbar=true,ShowMinIcethicknessLocations=true)
%
%%

arguments
    CtrlVar struct
    MUA     struct
    F       {mustBeA(F,{'struct','UaFields','numeric'})}
    Variable {mustBeA(Variable,{'string','numeric','logical'})}
    options.PlotGroundingLines  logical = true
    options.PlotCalvingFronts  logical = true
    options.CalvingFrontColor char = "b"
    options.GroundingLineColor char = "r"
    options.GetRidOfValuesDownStreamOfCalvingFronts=true;
    options.GetRidOfValuesDownStreamOfGroundingLines=false;
    options.PlotOverMesh=false;
    options.PlotUnderMesh=false;
    options.PlotMuaBoundary=true;
    options.ShowMinIcethicknessLocations=false;
    options.FigureTitle string="UaPlots";  % this is the figure title, not the plot title 
    options.CreateNewFigure logical = true ; 
    options.MeshColor char="k"
    options.logColorbar=false;
    options.Plot string = ""

    % options.ColorMap double=othercolor('YlGnBu6',1028)
    % options.ColorMap double=othercolor("Mlightterrain",1028)
    % options.ColorMap double=othercolor("Mdarkterrain",1028)
    % options.ColorMap double=othercolor("Mtemperaturemap",1028)
    % colormap(othercolor("Greys7",1028))
    % CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
    % CM=cmocean('ice',150) ; colormap(CM);
    
    options.ColorMap double=othercolor("YlGnBu8",1028)  % See othercolor.m for more options
end

%% Make F from old output files compatible

if ~isa(F,"UaFields")

    if  ~isfield(F,"LSF")
        F.LSF=[];
    end

    if  ~isfield(F,"x") || isempty(F.x)
        F.x=MUA.coordinates(:,1);
        F.y=MUA.coordinates(:,2);
    end

    if ~isfield(F,"time")
        F.time=[];
    end

    if ~isfield(F,"dt")
        F.dt=[];
    end

    if ~isfield(F,"GF")
        F.GF=[];
    elseif ~isfield(F.GF,"node")
        F.GF=[];
    end
else
    
   if isempty(F.x)
        F.x=MUA.coordinates(:,1);
        F.y=MUA.coordinates(:,2);
    end

end
%%

lg=[]; % this will be a handle to the legend (if created). 
%%

% if fig title has not been set, use by default the variable name
if options.FigureTitle=="UaPlots"
    if isstring(Variable)
       options.FigureTitle=Variable;
    elseif ~isempty(inputname(4))
        options.FigureTitle=inputname(4) ; 
    end

end


if options.CreateNewFigure
    fFig=FindOrCreateFigure(options.FigureTitle)  ; 
    
    clf(fFig)  ;
end

if isstring(Variable)

        options.Plot=Variable;

end

if islogical(Variable)
    Variable=double(Variable) ;
end

if isnumeric(Variable)

    % If the Variable is entered as a Nnodes x 2 array, then assume this is a velocity field and replace the (ub,vb) velocity in
    % F with this variable. This is an easy option to plot any velocity field over the mesh.
    [nV,mV]=size(Variable);
    if nV==MUA.Nnodes && mV==2
        F.ub=full(Variable(:,1));
        F.vb=full(Variable(:,2));
        if options.Plot==""
            options.Plot="-uv-" ;
        end

    else
        Variable=full(Variable);
    end
end

if isempty(F)
    F=UaFields;
end

if isstring(options.Plot)
% {"-eta-","eta int","etaint","-eta int-"}
    if contains(options.Plot,"int") || contains(options.Plot,"eta") || contains(options.Plot,"-e-") ...
            || contains(options.Plot,"tau")  || contains(options.Plot,"basal drag") 
        options.GetRidOfValuesDownStreamOfCalvingFronts=false;
    end

end

 

if options.GetRidOfValuesDownStreamOfCalvingFronts  && ~isempty(F.LSF)

    if isempty(F.LSFMask)
        F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    end

    F.ub(~F.LSFMask.NodesIn)=NaN;
    F.vb(~F.LSFMask.NodesIn)=NaN;

    if isnumeric(Variable)
        if numel(Variable)==MUA.Nnodes
            Variable(~F.LSFMask.NodesIn)=NaN;
        end
    end


end

if options.GetRidOfValuesDownStreamOfGroundingLines  && ~isempty(F.GF.node)  && options.Plot~="-strain rates-"

 

    F.ub(F.GF.node<0.5)=NaN;
    F.vb(F.GF.node<0.5)=NaN;

    if isnumeric(Variable)
        if numel(Variable)==MUA.Nnodes
            Variable(F.GF.node<0.5)=NaN;
        end
    end


end


xGL=nan ; yGL=nan ; xCF=nan ; yCF=nan ;

isModifyColormap=true;

if isModifyColormap
    colormap(options.ColorMap);
end

if options.PlotOverMesh
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    PlotMuaMesh(CtrlVar,MUA,nan,options.MeshColor) ;
    hold on

end


if isnumeric(Variable) && options.Plot ==""

    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Variable);
    title(cbar,inputname(4)) ;
    subtitle(sprintf("$t=%g \\quad  \\Delta t$=%g ",F.time,F.dt),Interpreter="latex")

else


    switch lower(options.Plot)

        case {"speed","-speed-"}

            speed=sqrt(F.ub.*F.ub+F.vb.*F.vb) ;

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\| \mathbf{v} \|$",Interpreter="latex")
            title(cbar,"$(\mathrm{m\,yr^{-1}})$",interpreter="latex")

        case {"log10speed","-log10speed-"}

            speed=sqrt(F.ub.*F.ub+F.vb.*F.vb) ;
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\log_{10}(\| \mathbf{v} \|)$",Interpreter="latex")
            title(cbar,["$\log_{10}(\| \mathbf{v} \|)$","(m/yr)"],Interpreter="latex")
            
            set(gca,'ColorScale','log')
            CM=cmocean('-ice',15) ; colormap(CM);

        case {"ubvb","-ubvb-","uv","-uv-"}

            CtrlVar.VelColorMap=jet(100) ;
            [cbar,~,CtrlVar]=QuiverColorGHG(F.x,F.y,F.ub,F.vb,CtrlVar) ;
            title(cbar,"(m/a)",Interpreter="latex")
            title("velocities",Interpreter="latex")
            subtitle(sprintf("$t=%g \\quad  \\Delta t$=%g ",F.time,F.dt),Interpreter="latex")


        case {"dhdt","-dhdt-","dh/dt","-dh/dt-"}


            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.dhdt);
            title(cbar,"(m/a)",Interpreter="latex")
            title(sprintf("$dh/dt$ at t=%g",CtrlVar.time),Interpreter="latex")
            title(cbar,"$(\mathrm{m\,yr^{-1}})$",interpreter="latex")

        case {"basal drag","taub"}


            [tbx,tby] = CalcBasalTraction(CtrlVar,[],MUA,F) ; 

            % [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,[],MUA,F) ;

            CtrlVar.VelColorMap=jet(100) ;
            [cbar,~,CtrlVar]=QuiverColorGHG(F.x,F.y,tbx,tby,CtrlVar) ;
            title(cbar,"(kPa)",Interpreter="latex")
            title(sprintf("basal drag vectors at t=%g",CtrlVar.time),Interpreter="latex")


        case "e node"  % effective strain rate


            [~,~,~,~,~,~,~,~,e,~]=CalcNodalStrainRatesAndStresses(CtrlVar,[],MUA,F) ;

            % e(e<0)=eps ; % the projection onto nodes does not preserve positive
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,e);
            title(cbar,"(1/a)",Interpreter="latex")
            title(sprintf("effective strain rates at t=%g",CtrlVar.time),Interpreter="latex")


        case {"-e-","e int","-e int-"}  % effective strain rate at integration points



            [~,~,~,~,~,~,~,e]=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n); % returns integration point values

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,e);
            title(cbar,"(1/a)",Interpreter="latex")
            title(sprintf("effective strain rates at integration points at t=%g",CtrlVar.time),Interpreter="latex")

        case "-strain rates-"

           [~,xint,yint,exx,eyy,exy]=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n); % returns integration point values

           if options.GetRidOfValuesDownStreamOfGroundingLines

               II=F.GF.ElementsDownstreamOfGroundingLines;
               exx(II,:)=0;
               eyy(II,:)=0;
               exy(II,:)=0;
    

           end

           
           scale=0.1 ; 
           LineWidth=1; 
           nStride=10;
           xint=xint(1:nStride:end,1);
           yint=yint(1:nStride:end,1);
           exx=exx(1:nStride:end,1);
           eyy=eyy(1:nStride:end,1);
           exy=exy(1:nStride:end,1);

           PlotTensor(xint/CtrlVar.PlotXYscale,yint/CtrlVar.PlotXYscale,exx,exy,eyy,scale,LineWidth)

        case "eta node"  % effective strain rate


            [~,~,~,~,~,~,~,~,~,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,[],MUA,F) ;

            % e(e<0)=eps ; % the projection onto nodes does not preserve positivy
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,eta);
            title(cbar,"(kPa yr)",Interpreter="latex")
            title(sprintf("effective viscosity eta at t=%g",CtrlVar.time),Interpreter="latex")


        case {"-eta-","eta int","etaint","-eta int-","eta"}  % effective strain rate at integration points



            etaInt=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n); % returns integration point values


            fFigHist=FindOrCreateFigure(options.FigureTitle+"Hist")  ; clf(fFigHist)  ; 
            histogram((log10(etaInt(:))),Normalization="probability") ; 
            hold on ; 
            xline(log10(CtrlVar.etaZero),'r',LineWidth=2)


            fFig=FindOrCreateFigure(options.FigureTitle)  ; clf(fFig)  ; 
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,log10(etaInt));
            title(cbar,"(kPa yr)",Interpreter="latex")
            title(sprintf("log10 of effective viscosity at integration points \n t=%g",CtrlVar.time),Interpreter="latex")



    

        case {"surface slope","-surface slope-"} 

            [dfdx,dfdy]=calcFEderivativesMUA(F.s,MUA,CtrlVar) ;
            slope=sqrt(dfdx.*dfdx+dfdy.*dfdy) ;
            slope=ProjectFintOntoNodes(MUA,slope) ;
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.rho.*F.h.*slope);
            title(cbar,"()",Interpreter="latex")
            title(sprintf("surface slope at t=%g",CtrlVar.time),Interpreter="latex")

        case {"log(ab)","-log(ab)-"}

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,-F.ab);
            set(gca,'ColorScale','log') 
            title(cbar,"$-a_b$ (m/yr)",Interpreter="latex")
            title("\textbf{Basal melt rates}",Interpreter="latex",FontSize=16)
            CM=cmocean('-thermal',20) ; CM(1,:)=[0.9 0.9 0.9]  ; colormap(CM);

        otherwise

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.(Variable));
            title(cbar,"$("+Variable+")$",Interpreter="latex")


    end
end

if options.logColorbar
    set(gca,'ColorScale','log')
end

hold on ;

if options.PlotUnderMesh
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    PlotMuaMesh(CtrlVar,MUA,nan,options.MeshColor) ;
    hold on

end


if options.PlotGroundingLines
    if isfield(F,"GF")
        [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],color=options.GroundingLineColor,DisplayName="Grounding lines");
    elseif isfield(F.GF,"node")
        [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,F.GF.node,[],[],[],color=options.GroundingLineColor,DisplayName="Grounding lines");
    end
end

if options.PlotCalvingFronts
    [xCF,yCF]=PlotCalvingFronts(CtrlVar,MUA,F,color=options.CalvingFrontColor);
end


if options.PlotMuaBoundary
    PlotMuaBoundary(CtrlVar,MUA,"b--",DisplayName="Mesh Boundary");
end



if isfield(CtrlVar,"PlotsXaxisLabel")
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")
    ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
end

if options.logColorbar
    set(gca,'ColorScale','log')
end

if options.ShowMinIcethicknessLocations

    iloc=F.h==CtrlVar.ThickMin ;
    nAtMinThick=numel(find(iloc));
    plot(F.x(iloc)/CtrlVar.PlotXYscale,F.y(iloc)/CtrlVar.PlotXYscale,LineStyle="none",Marker="o",MarkerFaceColor="m",MarkerEdgeColor="c",MarkerSize=3,DisplayName="at min thick")

    iloc=F.h<CtrlVar.ThickMin ;
    nLessThanMinThick=numel(find(iloc));
    plot(F.x(iloc)/CtrlVar.PlotXYscale,F.y(iloc)/CtrlVar.PlotXYscale,LineStyle="none",Marker="o",MarkerFaceColor="r",MarkerEdgeColor="c",MarkerSize=6,DisplayName="below min thick")

    [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.h,CtrlVar.ThickMin+eps(CtrlVar.ThickMin));
    plot(xc/CtrlVar.PlotXYscale,yc/CtrlVar.PlotXYscale,"g",LineWidth=2,DisplayName="Min thick")
     [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.h,2*CtrlVar.ThickMin);
    plot(xc/CtrlVar.PlotXYscale,yc/CtrlVar.PlotXYscale,"g",LineWidth=1,DisplayName="Twice min thick")
   
    
    title(sprintf("Nodes at/less than min thick (%i/%i)",nAtMinThick,nLessThanMinThick))
    lg=legend();

end


axis tight
subtitle(sprintf("$t=%g \\quad  \\Delta t$=%g ",F.time,F.dt),Interpreter="latex")

if ~nargout   % A trick to suppress any function output if no output requested. No need to suppress output using ;
    clearvars cbar
end




end