function cbar=UaPlots(CtrlVar,MUA,F,Variable,options)

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



if isnumeric(Variable)

    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Variable);

else


    switch lower(Variable)

        case {"speed","-speed-"}

            speed=sqrt(F.ub.*F.ub+F.vb.*F.vb) ;

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\| \mathbf{v} \|$",Interpreter="latex")

        case {"log10speed","-log10speed-"}

            speed=log10(sqrt(F.ub.*F.ub+F.vb.*F.vb)) ;
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\log_{10}(\| \mathbf{v} \|)$",Interpreter="latex")

        case {"ubvb","-ubvb-"}

            isModifyColormap=false;

%             if ~isempty(F.LSF)
%                 IC=F.LSF<0;
%                 F.ub(IC)=nan;
%                 F.vb(IC)=nan;
%             end
%             
            QuiverColorGHG(F.x,F.y,F.ub,F.vb,CtrlVar) ;

        otherwise

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.(Variable));


    end
end

hold on ;

if options.PlotGroundingLines
    PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],color=options.GroundingLineColor);
end

if options.PlotCalvingFronts
    PlotCalvingFronts(CtrlVar,MUA,F,color=options.CalvingFrontColor);
end

if isModifyColormap
colormap(options.ColorMap); 
end


end