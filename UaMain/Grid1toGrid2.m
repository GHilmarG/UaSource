function f2=Grid1toGrid2(DTgrid1,f1,x2,y2,CtrlVar,OutsideValue)
    
    % f2=Grid1toGrid2(DTgrid1,f1,x2,y2,CtrlVar,OutsideValue)
    % given values (f1) on nodal points on grid1 and the triangulation (DTgrid1) of grid 1
    % interpolate f1 onto (x2,y2)
    % if any of the x2 and y2 values contain NaN, f2 is returned with NaN at the corresponding locations
    % if OutsideValue is defined then any points outside of the convex hull are assigned that value
    % if OustideValue is not defined then these points are given the nearest-neighbor values
    %
    %
    %  Note: This is an old m-file created before MATLAB introduced the scatteredinterpoland object.
    %        Use MATLAB scatteredinterpolant instead, or consider using MapFbetweenMeshes.m
    % 
    %       
    % 

    %
    %%
    
    warning("Grid1toGrid2:ImTooOld","Old version, try not to use. Better options available.")

    if isempty(x2)
        f2=[];
        return
    end
    
    if nargin<6
        OutsideValue=NaN;
    end
    
    
    
    notnan=~isnan(x2) | ~isnan(y2) ;
    
    %     if any(isnan(x2)) || any(isnan(y2))
    %         fprintf(CtrlVar.fidlog,'nan in x y input fields to Grid1toGrid2  \n');
    %         save TestSave x2 y2
    %         f2=x2*0+nan;
    %         error('fdsa')
    %     end
    
    try
        F = TriScatteredInterp(DTgrid1,double(f1),'natural');
    catch err
        
        save TestSave
        rethrow(err);
        error(' error in TriScatteredIntep')
    end
    
    
    f2=x2*0;
    
    % I now interpolate where x2 an y2 are numbers
    f2(notnan) = F(x2(notnan),y2(notnan));
    
    % TriScatteredInterp returns NaN for all convex-set outside values
    % I now replace these with nearest neighbor values, or user defined value
    if any(isnan(f2))
        ind=isnan(f2)  & notnan ;
        
        if isnan(OutsideValue)
            nn=nearestNeighbor(DTgrid1,x2(ind),y2(ind));
            f2(ind)=f1(nn);
        else
            f2(ind)=OutsideValue;
        end
        
    end
    
    % finally I make sure that f2 is NaN where on input x2 and y2 where NaN
    f2(~notnan)=NaN;
    
    %     if any(isnan(f2))
    %         fprintf(CtrlVar.fidlog,'Grid1toGrid2 returns NaN  \n');
    %         save TestSave
    %         error('Error in Grid1Grid2: All variables were saved in file TestSave ')
    %     end
    
   
    
end

