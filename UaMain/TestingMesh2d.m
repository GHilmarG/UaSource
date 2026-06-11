

hdata = [];
options.output = false;
hdata.MeshSizeMax = 0.1;


node = [0 0; 1 0; 1 1; 0 1];              % Simple square example

[p,t] = mesh2d(node,[],hdata,options);    % Auto size fun only

figure
subplot(2,2,1)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Uniform mesh from MESH2D')

node = [0.3,0.3; 0.7,0.3; 0.7,0.7; 0.3,0.7];

in = inpoly(p,node);
ti = sum(in(t),2)>0;
[p,t] = refine(p,t,ti);

subplot(2,2,2)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Mesh refined in centre region using REFINE')

[p,t] = smoothmesh(p,t);

subplot(2,2,3)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Smoothed mesh using SMOOTHMESH')


node = [0.35,0.35; 0.65,0.35; 0.65,0.65; 0.35,0.65];

in = inpoly(p,node);
ti = sum(in(t),2)>0;
[p,t] = refine(p,t,ti);

subplot(2,2,4)
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Mesh refined in centre region using REFINE')


[p,t] = refine(p,t);[p,t] = smoothmesh(p,t);

figure
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','b');
axis equal off; hold on;
title('Mesh refined globally ')



