

%%

%load TestSave



  %%
  %
  %figure(2) ; hold on ; plot(GLele(:,[1 2])'/H0,GLele(:,[3 4])'/H0)
  %%
  %[GLx,GLy,GLxUpper,GLyUpper,GLxLower,GLyLower] = FindGL(DTxy,GF.node,500);
  
  %H0=1000; 
  %figure ; PlotFEmesh(coordinates/H0,connectivity,CtrlVar) ; axis tight ; 
  % hold on ; plot(GLx/H0,GLy/H0,'r');
  %%
GLele=GLgeometry(connectivity,coordinates,GF);
  H0=1000;
  hold on
  plot(GLele(:,[3 4])'/H0,GLele(:,[5 6])'/H0,'c','LineWidth',2)
  
%%
