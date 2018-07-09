
%%

CurDir=pwd;
goto_home_directory();

cd('Bamber\dtm')

load BamberVelDataGridded X Y Vx Vy Speed



%%




load BamberDataDelaunayTri DTxy x y vx vy speed


N=100;
x=x(1:N:end) ; y=y(1:N:end) ; speed=speed(1:N:end);

save TestSave x y speed

%%
load TestSave x y speed

%%
st=tpaps([x y]',speed');

%%

xmin=-2050 ; xmax=-1150 ; ymin=-750 ; ymax=100 ; dx=1 ; dy=1;

[X,Y] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

disp(' creating TriScatterInterpolant ')
F = TriScatteredInterp(DTxy, speed);  F.Method = 'natural'; Speed=F(X,Y);

cd(CurDir)

%%

