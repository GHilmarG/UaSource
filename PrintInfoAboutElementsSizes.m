



function [Emin,Emax,Emean,Emedian,Tlength]=PrintInfoAboutElementsSizes(CtrlVar,MUA,options)

%%
%
% Calculates equivalent element sizes
%
% Equivalent sizes can be defined in two ways:
%
% # as the leg of an isosceles right triangle with the same area as the triangular element, or
% # as the length of the side of a perfect square with the same area as the triangular element.
%
% Option 1) is the default, however after some thought it now feels to me that option 2) should be the default.
%
% These two different estimates is sqrt(2), with option (1) giving the numerically larger number.  
%
%
% Example:
%
%   [Emin,Emax,Emean,Emedian]=PrintInfoAboutElementsSizes(CtrlVar,MUA,print=true,plot=true,LengthMeasure="-side of a perfect square of equal area-");
%
%%

arguments
    CtrlVar struct
    MUA     struct
    options.print logical = true ;
    options.LengthMeasure string="-leg of an isosceles right triangle-" ;  % "-side of a perfect square of equal area-";
    options.plot logical = false ;
end


Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity);

if options.LengthMeasure=="-side of a perfect square of equal area-"
    Tlength=sqrt(Tarea)   ;   % this is the length of the sides of a perfect square of equal area.
else
    Tlength=sqrt(2*Tarea) ;   % Default option: this is the length of the legs of an isosceles right triangle.
end



Emax=max(Tlength);
Emean=mean(Tlength);
Emedian=median(Tlength);
Emin=min(Tlength);

if options.print

    fprintf(' #Elements: %-i, #Nod: %-i,  # Nodes=%-i. Elements have max, mean, median, and min sizes of %-g,  %-g,  %-g,  %-g,   respectively. \n',...
        MUA.Nele,MUA.nod,MUA.Nnodes,Emax,Emean,Emedian,Emin);
end

if options.plot

    cbar=UaPlots(CtrlVar,MUA,[],Tlength/CtrlVar.PlotXYscale) ;
    title(cbar,"$l$",Interpreter="latex")
    % Equivalent
    title("Equivalent element sizes",Interpreter="latex")

    if CtrlVar.PlotXYscale~=1

        StAdd="(/"+num2str(CtrlVar.PlotXYscale)+")" ;
        title(cbar,["$l$",StAdd],Interpreter="latex")
        title("Equivalent element sizes "+StAdd,Interpreter="latex")


    end


    FindOrCreateFigure("histogram of element sizes")
    histogram(Tlength/CtrlVar.PlotXYscale)

    title("Histogram of element sizes",Interpreter="latex")
    

    xlabel("Effective element size",Interpreter="latex")
    ylabel("Number of elements",Interpreter="latex")

    if CtrlVar.PlotXYscale~=1

        StAdd="(/"+num2str(CtrlVar.PlotXYscale)+")" ;
        title("Histogram of element sizes "+StAdd,Interpreter="latex")
        xlabel("Effective element size "+StAdd,Interpreter="latex")

    end

end


end
