function [Emin,Emax,Emean,Emedian]=PrintInfoAboutElementsSizes(CtrlVar,MUA)
    
    Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity); 
    Tlength=sqrt(2*Tarea) ;
    
    Emax=max(Tlength);
    Emean=mean(Tlength);
    Emedian=median(Tlength);
    Emin=min(Tlength);
    
    fprintf(' #Elements: %-i, #Nod: %-i,  # Nodes=%-i. Elements have max, mean, median, and min sizes of %-g,  %-g,  %-g,  %-g,   respectively. \n',...
        MUA.Nele,MUA.nod,MUA.Nnodes,Emax,Emean,Emedian,Emin);
    
    
end
