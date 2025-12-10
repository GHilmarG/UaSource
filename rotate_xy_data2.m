function [along,trans,angle]=rotate_xy_data2(x,y)

% fits best straigt line to xy coordinate data and rotates so that mean of 'trans'
% is zero and 'along' is the in-line displacement

x=x-x(1) ; y=y-y(1);
sol=[ones(length(x),1) x(:)]\[y(:)] ;
phi=atan(sol(2));

% test to see if the in-line displacement is positive 
% if not rotate the coordinate system by 180 degrees
along=x(end)*cos(phi)+y(end)*sin(phi);
if along < 0
  phi=phi+pi; 
end

along=x*cos(phi)+y*sin(phi);
trans=-x*sin(phi)+y*cos(phi);



angle=phi*180/pi;

end
