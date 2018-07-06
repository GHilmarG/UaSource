function SetRegionalPlotAxis(Region)


%%
% SetRegionalPlotAxis(Region)
%
% Region={pigiceshelf|thwaitesiceshelf|pig-thwaites|filchner-ronne
%
%

switch lower(Region)
    
    case 'pigiceshelf'
        
        axis([-1680 -1540 -380 -240])
        
    case 'thwaitesiceshelf'
        
        axis([-1620 -1500 -520 -400])
        
    case 'pig-twg-shelves'
        
        axis([-1700 -1500 -510 -250])
        
    case {'pig-thwaites','pt','pig-twg'}
        
        axis([-1720 -1350 -520 -210])
    case{'filchner-ronne','fr','rf'}
        axis([-1600 -400 100 1100])
    case'ross'
        axis([-600 400 -1400 -400])
end


end


