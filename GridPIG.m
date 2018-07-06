

cd ../../TimmermannDataSet

M=dlmread('RTopo104_gl.asc');
gl.lon=M(:,1); gl.lat=M(:,2);
[gl.x,gl.y]=geog_to_pol_wgs84_71S(gl.lat,gl.lon);

%figure ; plot(gl.x,gl.y,'.')


M=dlmread('RTopo104_coast.asc');
coast.lon=M(:,1); coast.lat=M(:,2);

ind=coast.lat<=-60;  coast.lon=coast.lon(ind) ; coast.lat=coast.lat(ind);
[coast.x,coast.y]=geog_to_pol_wgs84_71S(coast.lat,coast.lon);

%figure ; plot(coast.x,coast.y,'r.') ; axis equal


fprintf(' loading topo data ')
load TGridHeightPIG xequal yequal sgrid x y s 
load TGridDraftPIG xequal yequal bgrid x y b 
load TGridBathyPIG xequal yequal Bgrid x y B 
fprintf(' topo data loaded \n')

cd ../ssa/FEicestream2d


hgrid=sgrid-bgrid;


figure
contourf(xequal,yequal,hgrid,50)
axis equal
colorbar
title('h')



%exctract zero thickness contour
hgrid(1,:)=0 ; hgrid(end,:)=0 ; hgrid(:,1)=0 ; hgrid(:,end)=0;
[xc,yc]=getcon(xequal,yequal,hgrid,0);
figure ; axis equal
plot(xc,yc) ; axis equal
ind=find(isnan(xc)); ind=[1;ind] ;
hold on
col=['b','r','g','c','k'] ;
for I=1:length(ind)-1
	plot(xc(ind(I):ind(I+1)),yc(ind(I):ind(I+1)),'color',col(I))
end

% from visual inspection I see I need the gree curve

hnillx=xc(ind(3)+1:ind(4)-1); hnilly=yc(ind(3)+1:ind(4)-1);
hold on;
plot(hnillx,hnilly,'k')

N=10;
xy=[hnillx(1:N:end) hnilly(1:N:end)];

%%
MaxEleSize=100e3;
DesiredEleSize=100e3;

hfunOpt=DesiredEleSize;
hdata.fun = @hfunPIG;
hdata.args  = {hfunOpt} ;
hdata.MeshSizeMax  = MaxEleSize; 

[p,t] = mesh2d(xy,[],hdata);




%%



figure
contourf(xequal,yequal,sgrid,50)
axis equal
colorbar
title('s')

%%
figure
contourf(xequal,yequal,Bgrid,50)
axis equal
colorbar
title('B')
%%

figure
plot(coast.x,coast.y,'b.') ; hold on
plot(gl.x,gl.y,'r.')

%%




%%



