function dL=TidalMigrationDistance(rho,rhow,dhdx,dsdx,dBdx,dS)

% calculates the GL migration distance for a given shift in sea level 
% or a given tidal amplitude
%
% Two cases are considered:
% 1) dhdx is fixed across the grounding line
% 2) dsdx is fixed across the grounding line
%
% To get case (1) leave dsdx empty, for case (2) leave dhdx empty
%
% On return dL is the shift in grounding line position for 
% low (negative) and high (positive) tide.  Give dS (tidal amplitude) as a positive number
%

dS=abs(dS);

if isempty(dsdx)
    
    dL.pos=rhow*dS/(rho*dhdx+rhow*dBdx);
    dL.neg=-dL.pos;
    
elseif isempty(dhdx)
    
    %dhdx=dsdx/(1-rho/rhow);
    %dL.neg=-rhow*dS/(rho*dhdx+rhow*dBdx);
    
    %dhdx=dsdx-dBdx;
    %dL.pos=rhow*dS/(rho*dhdx+rhow*dBdx);
    
    dL.neg=-dS/( ((rho/rhow)/(1-rho/rhow)) *dsdx+dBdx);
    dL.pos=-dL.neg/(1-rho/rhow);
    
    
else
    
    fprintf('either dhdx or dsdx must be empty \n')
    dL.neg=NaN;
    dL.pos=NaN;
    
end

end
