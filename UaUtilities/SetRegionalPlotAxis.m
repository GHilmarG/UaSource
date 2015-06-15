function SetRegionalPlotAxis(Region)

switch lower(Region)
    
    case 'pigiceshelf'
        
        axis([-1680 -1540 -380 -240])
        
    case 'thwaitesiceshelf'
        
        axis([-1620 -1500 -520 -400])
       
    case 'pig-thwaites' 
        
        axis([-1720 -1350 -520 -210])
        
end


end


