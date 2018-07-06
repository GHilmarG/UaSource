function goto_home_directory()
	
	% tries to find my home directory
    
	% is `ghg' or `hilmar' a part of current path?
	wd=pwd;
    indsGHG=strfind(upper(wd),'SETTINGS\GHG');
	indGHG=strfind(upper(wd),'GHG');
	indHilmar=strfind(upper(wd),'HILMAR');
    indsHilmar=strfind(upper(wd),'SETTINGS\HILMAR');
	
   
    if ~isempty(indGHG) && isempty(indsGHG)
        hd=wd(1:indGHG+2);
        cd(hd);
    elseif ~isempty(indHilmar) && isempty(indsHilmar)
        hd=wd(1:indHilmar+5);
        cd(hd);
    else
       
        try
            cd /media/Master2011/GHG
        catch
            try
                cd('E:\GHG')
            catch
                fprintf(' can not find home directory')
            end
        end
        
    end
    
end
%%










