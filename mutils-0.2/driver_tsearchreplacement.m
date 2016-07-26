x=[0 1 0 1 ]; y=[0 0 1 1];

TRI=delaunay(x,y);

xi=[-0.25; 0.75] ; yi=[0.25; 0.75] ;

T1=tsearch(x,y,TRI,xi,yi)

T2=tsearchReplacement(x,y,TRI,xi,yi)

