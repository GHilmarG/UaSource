function cbar=UaPlots(CtrlVar,MUA,F,Variable,options)

%%
%
% cbar=UaPlots(CtrlVar,MUA,F,Variable,options)
%
% Simple plot utility to plot variables and calving fronts and grounding lines as well.
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
%%

arguments
    CtrlVar struct
    MUA     struct
    F       {mustBeA(F,{'struct','UaFields'})}
    Variable   {string, double}
    options.PlotGroundingLines  logical = true
    options.PlotCalvingFronts  logical = true
    options.CalvingFrontColor char = "k"
    options.GroundingLineColor char = "r"
    options.GetRidOfValuesDownStreamOfCalvingFronts=true;
    options.PlotOverMesh=false;
    options.PlotUnderMesh=false;


    % options.ColorMap double=othercolor('YlGnBu6',1028)
    % options.ColorMap double=othercolor("Mlightterrain",1028)
    % options.ColorMap double=othercolor("Mdarkterrain",1028)
    % options.ColorMap double=othercolor("Mtemperaturemap",1028)
    options.ColorMap double=othercolor("YlGnBu8",1028)  % See othercolor.m for more options
end


if islogical(Variable)
    Variable=double(Variable) ;
end

if options.GetRidOfValuesDownStreamOfCalvingFronts  && ~isempty(F.LSF)

    if isempty(F.LSFMask)
        F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    end

    F.ub(~F.LSFMask.NodesIn)=NaN;
    F.vb(~F.LSFMask.NodesIn)=NaN;

end

isModifyColormap=true;

if isModifyColormap
    colormap(options.ColorMap);
end

if options.PlotOverMesh
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    PlotMuaMesh(CtrlVar,MUA) ;
    hold on 

end


if isnumeric(Variable)

    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Variable);
    title(cbar,inputname(4)) ; 

else


    switch lower(Variable)

        case {"speed","-speed-"}

            speed=sqrt(F.ub.*F.ub+F.vb.*F.vb) ;

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\| \mathbf{v} \|$",Interpreter="latex")
            title(cbar,"$(\mathrm{m\,yr^{-1}})$",interpreter="latex")

        case {"log10speed","-log10speed-"}

            speed=log10(sqrt(F.ub.*F.ub+F.vb.*F.vb)) ;
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\log_{10}(\| \mathbf{v} \|)$",Interpreter="latex")
            title(cbar,"$\log_{10}(m/a)$",Interpreter="latex")
      

        case {"ubvb","-ubvb-"}

            CtrlVar.VelColorMap=jet(100) ; 
            cbar=QuiverColorGHG(F.x,F.y,F.ub,F.vb,CtrlVar) ;
            title(cbar,"(m/a)")
            title(sprintf("velocities at t=%f",CtrlVar.time),Interpreter="latex")

        otherwise

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.(Variable));
            title(cbar,"$(m)$",Interpreter="latex")


    end
end

hold on ;

if options.PlotUnderMesh
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    PlotMuaMesh(CtrlVar,MUA,[],"w") ;
    hold on 

end


if options.PlotGroundingLines
    PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],color=options.GroundingLineColor);
end

if options.PlotCalvingFronts
    PlotCalvingFronts(CtrlVar,MUA,F,color=options.CalvingFrontColor);
end



% Just guessing that this might be the most common case, the user can easily change afterwards anyhow.
xlabel("xps (km)",Interpreter="latex")
ylabel("yps (km)",Interpreter="latex")


end