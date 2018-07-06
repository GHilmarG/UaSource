function [b] = PIGfbed(arguments)
    
    %  Gives glacier lower surface for PIG for any (x,y) values
    %
    
    
    xequal=arguments{1};
    yequal=arguments{2};
    bgrid=arguments{3};
    coordinates=arguments{4};
    
    b = interp2(xequal,yequal,bgrid,coordinates(:,1),coordinates(:,2),'linear');
    
    
end
