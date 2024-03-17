function [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,varargin)

%
% Plots calving fronts, just a simple wrapper around PlotGroundingLines
%
%
% To suppress plotting set CtrlVar.PlotGLs=false
%
%
% F is the Ua field variable, but it can also be given as LSF or any other nodal variable
%
% Example:
%
% figure ; [xC,yC]=PlotCalvingFronts(CtrlVar,MUA,F,color="r");
%
% figure ; [xC,yC]=PlotCalvingFronts(CtrlVar,MUA,LSF,color="r");
%
%
%  figure ; [xC,yC]=PlotCalvingFronts(CtrlVar,MUA,LSF,"k--");
%
%
%  figure ; PlotCalvingFronts([],"ITS-LIVE",[],"r");
%
% Also consider using:
%
%   [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,Field,Value,varargin)
%
%
% Note: If called with not arguments, or with the second argument equal to "ITS-LIVE" , the outlines of the Antarctic Ice
% Sheet are plotted based on the ITS-LIVE120  ocean mask
%
%%



if nargin==0

    CtrlVar=[];
    MUA="ITS-LIVE" ;

end

if isempty(CtrlVar)
    CtrlVar(1).PlotXYscale=1000;
end


if isstring(MUA)

    GroundingLineDataSet=MUA ;

    switch GroundingLineDataSet

        case "ITS-LIVE"

            DataFile="ITS-LIVE-ANT-G0120-0000-AntarticIceSheetBoundary-nStride5.mat";

            if isfile(DataFile)
                load(DataFile,"AISBoundaryITS") ;
                xc=AISBoundaryITS(:,1);
                yc=AISBoundaryITS(:,2);
            else
                fprintf("Can not plot calving fronts based on ITS-LIVE data because the data file %f \n",DataFile)
                fprintf("is not found.\n")
                xc=[] ; yc=[]; 
            end



            otherwise

            error("case not found")

    end

    tt=axis;
    if isempty(varargin)
        plot(xc/CtrlVar.PlotXYscale,yc/CtrlVar.PlotXYscale,'b') ;
    else
        plot(xc/CtrlVar.PlotXYscale,yc/CtrlVar.PlotXYscale,varargin{:}) ;
    end
    ax=gca; ax.DataAspectRatio=[1 1 1];

    if ~isequal(tt,[0 1 0 1])
        axis(tt)
    end

    return

end



if  isnumeric(F) && numel(F)==MUA.Nnodes
    LSF=F;
else
    LSF=F.LSF ;
end

if isempty(LSF)
    xc=[] ; yc=[] ;
    return
end


CtrlVar.LineUpGLs=true ;

GF.node=LSF ;

CtrlVar.GLthreshold=0;
[xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],varargin{:}) ;





end