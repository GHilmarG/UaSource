function [xps,yps,u,v]=GetMeasuresVelocities(DataSet,Region)

%
%  DataSet and Region are optional inputs
%  
%  DataSet={[],'990m','450m'}
%  The boundary of the data set can be defined ast
%   Region.xmin <= x <= Region.xmax 
%   Region.ymin <= y <= Region.ymax
%
%
%

persistent  x y vx vy 

if nargin==0 || isempty(DataSet)
    DataSet='990m';
end

if isempty(x)
    
    AntarcticGlobalDataSets=getenv('AntarcticGlobalDataSets');
    
    if isempty(AntarcticGlobalDataSets)
        error('The environmental variable AntarcticDataSets not defined' )
    end
    
    switch DataSet
        
        case '990m'
            fprintf(' loading MEASURES 990m velocity data ')
            locdir=pwd;
            cd(AntarcticGlobalDataSets)
            cd MEASURES/990m
            load MEASURES_990m_x_y_vx_vy_err.mat x y vx vy err
            fprintf(' done \n ')
            cd(locdir)
        case '450m'
            fprintf(' loading MEASURES 450m velocity data ')
            locdir=pwd;
            cd(AntarcticGlobalDataSets)
            cd MEASURES/450m
            load MEASURES_450m_x_y_vx_vy_err.mat x y vx vy err
            fprintf(' done \n ')
            cd(locdir)
        otherwise
            error(' which case?')
            
    end
    
    vx=rot90(vx) ; vy=rot90(vy);
    % in Eric's data set, no meas are replaced by zero (!)
    vx(vx==0)=NaN; vy(vy==0)=NaN;  % here put no values to NaN
end

if nargin>1 && ~isempty(Region)
    if numel(Region)==4
        Ix=x>= Region(1)  & x <= Region(2);
        Iy=y>= Region(3)  & y <= Region(4);
    else
        Ix=x>= Region.xmin  & x <= Region.xmax;
        Iy=y>= Region.ymin  & y <= Region.ymax;
    end
    
    xps=x(Ix) ; yps=y(Iy) ; u=vx(Iy,Ix) ; v=vy(Iy,Ix);
else
    xps=x ; yps=y ; u=vx ; v=vy;
end


end

