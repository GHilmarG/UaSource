function [uMeas,vMeas,Err,Mask]=EricVelocities(CtrlVar,coordinates,DataSet)
    
    %
    % Mask is true where returned values are zero, ie where there is no underlying data available
    %
    
    persistent  x y vx vy Fvx Fvy Ferr
    
    if nargin<3
        DataSet='990m';
    end
    
    if isempty(x)
        
        AntarcticGlobalDataSets=getenv('AntarcticGlobalDataSets');
%         
%         if isempty(AntarcticGlobalDataSets)
%             error('The environmental variable AntarcticDataSets not defined' )
%         end
        
        switch DataSet
            case 'Eric'  % possible the same as 990m, a bit older version
                fprintf(' loading Eric velocity data ')
                locdir=pwd;
                cd(AntarcticGlobalDataSets)
                cd Antarctic' Ice Velocity Rignot'/
                load VelRignot_x_y_vx_vy_err.mat x y vx vy err
                fprintf(' done \n ')
                
                
            case '990m'
                fprintf(' loading MEASURES 990m velocity data ')
                locdir=pwd;
                
                if ~isempty(AntarcticGlobalDataSets)
                    cd(AntarcticGlobalDataSets)
                    cd MEASURES/990m
                end
                
                load MEASURES_990m_x_y_vx_vy_err.mat x y vx vy err
                fprintf(' done \n ')
                
            case '450m'
                fprintf(' loading MEASURES 450m velocity data ')
                locdir=pwd;
                cd(AntarcticGlobalDataSets)
                cd MEASURES/450m
                load MEASURES_450m_x_y_vx_vy_err.mat x y vx vy err
                fprintf(' done \n ')
                
            otherwise
                error(' which case?')
                
        end
        cd(locdir)
        
        % in Eric's data set, no meas are replaced by zero (!)
        vx(vx==0)=NaN; vy(vy==0)=NaN;  % here put no values to NaN
        
        err=single(err); % this is byte signed in the original data set
        
        [X,Y]= ndgrid(x,y);
        Fvx = griddedInterpolant(X,Y,flipud(rot90(vx,2)), 'linear');
        Fvy = griddedInterpolant(X,Y,flipud(rot90(vy,2)), 'linear');
        Ferr = griddedInterpolant(X,Y,flipud(rot90(err,2)), 'linear');
        
    end
    
    %speed=sqrt(vx.*vx+vy.*vy);
    %N=100 ; figure ; contourf(x(1:N:end)/1000,y(1:N:end)/1000,(rot90(log10(speed(1:N:end,1:N:end)),1)),'LineStyle','none') ; axis equal
    %axis equal ; colorbar
    
  
        
    % Because I put NaN in where there is no data, this will give NaN at location not surounded by four data points
    % I can then afterwards find the NaN and put in some data with very high errors
    uMeas=Fvx(coordinates(:,1),coordinates(:,2));
    vMeas=Fvy(coordinates(:,1),coordinates(:,2));
    Err=Ferr(coordinates(:,1),coordinates(:,2));
    uMeas=double(uMeas) ; vMeas=double(vMeas); Err=double(Err);
    %Err=zeros(length(coordinates),1)+1;
    
    
    % I need values everywhere, so I now set vel to zero where no measurments are available, 
    % and the error to a very large value
    Mask=isnan(uMeas) | isnan(vMeas);
    uMeas(Mask)=0 ; vMeas(Mask)=0;
    Err(Mask)=1e10;
    
   
    
    
end

