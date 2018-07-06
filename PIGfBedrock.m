function [B] = PIGfBedrock(arguments)
    
    %  Gives Bedrock for PIG for any (x,y) values
    %
    
    
    xequal=arguments{1};
    yequal=arguments{2};
    Bgrid=arguments{3};
    coordinates=arguments{4};
    
    B = interp2(xequal,yequal,Bgrid,coordinates(:,1),coordinates(:,2),'linear');
    
    
end
