function CompareResultsWithPreviouslyObtainedResults(Experiment,coordinates,connectivity,u,v,s,b,S,B,time,dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar)

[u0,v0,w0,s0,b0,S0,B0,coordinates0,connectivity0,DTxy0,TRIxy0] = loadCompareFile();

 

x=coordinates(:,1); y=coordinates(:,2); h=s-b;

% interpolate comparison dat on numerical grid

s0=Grid1toGrid2(DTxy0,s0,x,y);
u0=Grid1toGrid2(DTxy0,u0,x,y);
v0=Grid1toGrid2(DTxy0,v0,x,y);
w0=Grid1toGrid2(DTxy0,w0,x,y);

% calculate numerical vertical velocities
%


[etaInt,xint,yint,exx,eyy,exy,Eint,e]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
[w,wSurfInt,wBedInt]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,u,v,exx,eyy,xint,yint,coordinates,connectivity,nip,CtrlVar);



figure;
subplot(3,4,1) ; trimesh(TRIxy,x,y,s); title('s') ; xlabel('x') ; ylabel('y')
subplot(3,4,5); trimesh(TRIxy,x,y,s0); title('s (comparison)') ; xlabel('x') ; ylabel('y')
subplot(3,4,9); trimesh(TRIxy,x,y,s-s0); title('s difference ') ; xlabel('x') ; ylabel('y')

subplot(3,4,2); trimesh(TRIxy,x,y,u); title('u  ') ; xlabel('x') ; ylabel('y')
subplot(3,4,6); trimesh(TRIxy,x,y,u0); title('u (comparison)') ; xlabel('x') ; ylabel('y')
subplot(3,4,10); trimesh(TRIxy,x,y,u-u0); title('u difference') ; xlabel('x') ; ylabel('y')

subplot(3,4,3); trimesh(TRIxy,x,y,v); title('v  ') ; xlabel('x') ; ylabel('y')
subplot(3,4,7); trimesh(TRIxy,x,y,v0); title('v (comparison)') ; xlabel('x') ; ylabel('y')
subplot(3,4,11); trimesh(TRIxy,x,y,v-v0); title('v difference') ; xlabel('x') ; ylabel('y')

subplot(3,4,3); trimesh(TRIxy,x,y,w); title('w  ') ; xlabel('x') ; ylabel('y')
subplot(3,4,7); trimesh(TRIxy,x,y,w0); title('w (comparison)') ; xlabel('x') ; ylabel('y')
subplot(3,4,11); trimesh(TRIxy,x,y,w-w0); title('w difference') ; xlabel('x') ; ylabel('y')




u=u-mean(u) ; u0=u0-mean(u0);
v=v-mean(v) ; v0=v0-mean(v0);
w=w-mean(w) ; w0=w0-mean(w0);

diffs=norm(s-s0)/sqrt(length(s));
diffu=norm(u-u0)/sqrt(length(u));
diffv=norm(v-v0)/sqrt(length(v));
diffw=norm(w-w0)/sqrt(length(w0));
%             diffw=norm(wSurfInt(:)-w0(:))/sqrt(length(w0(:)));

disp(['         rms s: ',num2str(diffs)])
disp(['         rms u: ',num2str(diffu)])
disp(['         rms v: ',num2str(diffv)])
disp(['  rms w (node): ',num2str(diffw)])
%             disp([' rms w: ',num2str(diffw)])

[Nele,nod]=size(connectivity);

Nnodes=length(s);
fid = fopen('CompareNumericalNumerical.txt', 'a+');
fprintf(fid, '%i \t %i \t %g \t %g \t %g \t %g \t %i \t %i \n', Nele,nod,diffs,diffu,diffv,diffw,nip,Nnodes);
fclose(fid);


% view the contents
end

%   Gauss b pert     100x100km alpha=0.001; hmean=1000; ampl_b=100; sigma_bx=10000 ; sigma_by=10000;
%    Nele        nod       u rms     v rms    w rms
%    2048         6
%
%

