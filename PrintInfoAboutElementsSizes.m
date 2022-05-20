function [Emin,Emax,Emean,Emedian]=PrintInfoAboutElementsSizes(CtrlVar,MUA,options)


arguments
    CtrlVar struct
    MUA     struct
    options.print logical = true ;

end


Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity);
Tlength=sqrt(2*Tarea) ;

Emax=max(Tlength);
Emean=mean(Tlength);
Emedian=median(Tlength);
Emin=min(Tlength);

if options.print

    fprintf(' #Elements: %-i, #Nod: %-i,  # Nodes=%-i. Elements have max, mean, median, and min sizes of %-g,  %-g,  %-g,  %-g,   respectively. \n',...
        MUA.Nele,MUA.nod,MUA.Nnodes,Emax,Emean,Emedian,Emin);
end

end
