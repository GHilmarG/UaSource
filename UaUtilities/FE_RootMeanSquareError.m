function FERMSE=FE_RootMeanSquareError(f,g,M,h)
    
%
% [FERMSE,M]=FE_RootMeanSquareError(f,g,M,Normalize)
%
% Calculate the root-mean-square deviation/error between the nodal variables f and g.
%
% This is based on the inner-produced induced norm L2
%
% If h is given as an input the results is normalized by the norm of h. 
% 
%
% Example:
%
% To calculate FE-RMSE between f and g and normalize by f:
%
%   [FERMSE,M]=FE_RootMeanSquareError(f,g,MUA.M,f)
%
% 
% To calculate FE-RMSE between f and g and normalize by area.
%
%   [FERMSE,M]=FE_RootMeanSquareError(f,g,MUA.M,f*0+1)
%

    
    FERMSE=sqrt(FE_inner_product(f-g,f-g,M)) ;
    
    if nargin==4
        FERMSE=FERMSE/sqrt(FE_inner_product(h,h,M)) ;
    end
    
    
end

