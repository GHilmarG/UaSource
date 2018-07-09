function [xequal,yequal,sgrid,bgrid,Bgrid,MeshBoundaryCoordinates] = PIGdata(N)
	
	InfoLevel=10;
	
    	xmin=-2000e3 ; xmax=-1250e3 ; ymin=-500e3 ; ymax=0e3;
    
	%
	% returns Surface (s), Glacer bed  (b), and Ocean floor (B) on an rectangular grid for PIG
	%
	
	%N=50;  % only every Nth point of the zero thickness contour is used to define the boundary
	
	
	cd ../../TimmermannDataSet
	M=dlmread('RTopo104_coast.asc');
	coast.lon=M(:,1); coast.lat=M(:,2);
	
	ind=coast.lat<=-60;  coast.lon=coast.lon(ind) ; coast.lat=coast.lat(ind);
	[coast.x,coast.y]=geog_to_pol_wgs84_71S(coast.lat,coast.lon);
	
	
	fprintf(' loading topo data ')
	load TGridHeightPIG xequal yequal sgrid x y s
	load TGridDraftPIG xequal yequal bgrid x y b
	load TGridBathyPIG xequal yequal Bgrid x y B
	fprintf(' topo data loaded \n')
	
	
% 	% for some reason the boundaries have NaN
% 	
% 	xequal=xequal(2:end-1); yequal=yequal(2:end-1);
% 	
% 	sgrid=sgrid(2:end-1,2:end-1);
% 	bgrid=bgrid(2:end-1,2:end-1);
% 	Bgrid=Bgrid(2:end-1,2:end-1);
% 	
	
	cd ../ssa/FEicestream2d
	
	
	hgrid=sgrid-bgrid;
	
	%exctract zero thickness contour
	hgrid(1,:)=0 ; hgrid(end,:)=0 ; hgrid(:,1)=0 ; hgrid(:,end)=0;
	[xc,yc]=getcon(xequal,yequal,hgrid,0);
	
	ind=find(isnan(xc)); ind=[1;ind] ;
	
	if InfoLevel>1
		col=['b','r','g','c','k'] ;
		for I=1:length(ind)-1
			plot(xc(ind(I):ind(I+1)),yc(ind(I):ind(I+1)),'color',col(I))
		end
	end
	
	% from visual inspection I see I need the green curve
	hnillx=xc(ind(3)+1:ind(4)-1); hnilly=yc(ind(3)+1:ind(4)-1);
	
	
	% I must exclude points outside of the desired model domain, which in this case is a simple square from
	
	

	
	
	ind=(hnillx >= xmin ) & (hnillx <= xmax ) & (hnilly>= ymin ) & (hnilly<= ymax) ;
	
	hnillx=hnillx(ind) ; hnilly=hnilly(ind);
	
	
	Iterations=10;
	hnillx=SimpleMeanSmooth(hnillx,Iterations);
	hnilly=SimpleMeanSmooth(hnilly,Iterations);
	
% 	for I=1:50
% 		tempx=((hnillx(1:end-2)+hnillx(3:end))/2+hnillx(2:end-1))/2;
% 		tempy=((hnilly(1:end-2)+hnilly(3:end))/2+hnilly(2:end-1))/2;
% 		hnillx(2:end-1)=tempx ; hnilly(2:end-1)=tempy ;
% 	end
% 	
	MeshBoundaryCoordinates=[hnillx(1:N:end) hnilly(1:N:end)];
	
	
	% Smooth boundary nodes
	
	
	
	
	[n1,n2]=size(MeshBoundaryCoordinates);
	MeshBoundaryCoordinates(n1,1)=MeshBoundaryCoordinates(n1,1) ; MeshBoundaryCoordinates(n1,2)=ymax ;
	MeshBoundaryCoordinates(n1+1,1)=xmax ; MeshBoundaryCoordinates(n1+1,2)=ymax ;
	MeshBoundaryCoordinates(n1+2,1)=xmax ; MeshBoundaryCoordinates(n1+2,2)=ymin ;
	MeshBoundaryCoordinates(1,1)=MeshBoundaryCoordinates(1,1) ; MeshBoundaryCoordinates(1,2)=ymin ;
	
	
	
	if InfoLevel>1
		figure ; plot(MeshBoundaryCoordinates(:,1)/1000,MeshBoundaryCoordinates(:,2)/1000,'o-')
	end
	
	
	
end