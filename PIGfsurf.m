function [s] = PIGfsurf(arguments)
    
    %  Gives Bedrock for PIG for any (x,y) values
    %
    
    
    xequal=arguments{1};
    yequal=arguments{2};
    sgrid=arguments{3};
    coordinates=arguments{4};
    
    s = interp2(xequal,yequal,sgrid,coordinates(:,1),coordinates(:,2),'linear');
    
    
end
