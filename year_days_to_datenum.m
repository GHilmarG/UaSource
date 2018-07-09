function daynr=year_days_to_datenum(year,dayofyear,hour,min,sec)

% calculates matlap daynr from year and dayofyear
% for example: 
% datestr(year_days_to_datenum(2007,1))='01-Jan-2007'
% datestr(year_days_to_datenum(2007,0))='31-Dec-2006'



if nargin ==2 
  daynr=datenum(year,0,0)+dayofyear;
elseif nargin==5 
  daynr=datenum(year,0,0)+dayofyear+hour/24+min/24/60+sec/24/60/60;
else
  disp(' incorrect nr of arguments ') ; 
end
  
return
