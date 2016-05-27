function PlaceLabel(pt,str)

% simple function to place lables on plots (stole this from the internet somewhere)

x = pt(1);
y = pt(2);
h = line(x,y);
h.Marker = '.';
h = text(x,y,str);
h.HorizontalAlignment = 'center';
h.VerticalAlignment = 'bottom';


end
