function PrintInfoAboutElementsSizes(CtrlVar,MUA)
    
    Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity); 
    Tlength=sqrt(2*Tarea) ;
    
    fprintf(' #Elements: %-i, #Nod: %-i,  # Nodes=%-i. Elements have max, mean, median, and min sizes of %-g,  %-g,  %-g,  %-g,   respectively. \n',...
        MUA.Nele,MUA.nod,MUA.Nnodes,max(Tlength),mean(Tlength),median(Tlength),min(Tlength));
    
    
end
