
function GLQ=GroundingLineCurvePlots2(UserVar,CtrlVar,MUA,F,options)


%% GroundingLineCurvePlots2("0000050-Nodes83632-Ele166223-Tri3-kH10000-FT-P-TWIS-MR4-SM-5km-Alim-Ca1-Cs100000-Aa1-As100000-")
% GroundingLineCurvePlots2("D:\Runs\Calving\PIG-TWG\ResultsFiles\0000050-Nodes83632-Ele166223-Tri3-kH10000-FT-P-TWIS-MR4-SM-5km-Alim-Ca1-Cs100000-Aa1-As100000-")



arguments
    UserVar struct
    CtrlVar struct
    MUA     struct
    F       struct
    options.xgl   (1,:)  {mustBeNumeric}=NaN
    options.ygl   (1,:)  {mustBeNumeric}=NaN
    options.smoothing  (1,1) {mustBeNumeric}=1 % i.e. by default no smoothing, this value is between 0 and 1
    options.ds  (1,1) {mustBeNumeric}=50e3 % distance between points along GL for which quantities are calculated
    options.Plot logical=false 

end


MUA=UpdateMUA(CtrlVar,MUA);  % Maybe MUA needs updating,


%%  1) Get stresses
[txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,UserVar,MUA,F) ;
F.txx=txx ;
F.txy=txy ;
F.tyy=tyy ;

%% 2) now get the grounding line or grounding-line segment of interest, or possibly this was given on input
% and calculate quantities

colormap('HSV')


if isnan(options.xgl)

    if options.Plot
        CtrlVar.PlotGLs=1;
        FigGL=FindOrCreateFigure("Grounding lines") ;
    else
        CtrlVar.PlotGLs=0;
    end

    [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,F.GF) ;
    % here select just one grounding line
    I=find(isnan(xGL)) ;
    xgl=xGL(1:I(1));  ygl=yGL(1:I(1));


else

    xgl=options.xgl;
    ygl=options.ygl;
end



% Get points along the grounding-line segment at equally spaced intervals
% It is quite possible that some smoothing of the grounding line will here be requried
% The GL calculated tends to be locally quite sensitive to the mesh resolution and
% the shapes of eleents.


%% 3) Now sample grounding-line at equal distance ds, and possibly smooth the grounding-lines as well
% smoothing=1e-9 ;ds=2e3 ;

[xglc,yglc,nxGL,nyGL,tvector] = Smooth2dPos(xgl,ygl,options.smoothing,options.ds) ;



%%  4) Here the GL quantities are calculated
GLQ=GroundingLineQuantities(F,xglc,yglc,nxGL,nyGL) ;
GLQ.xglc=xglc;
GLQ.yglc=yglc;
GLQ.nxGL=nxGL;
GLQ.nyGL=nyGL;


% Return? 

if ~options.Plot
    return
end


%--------------------------------------------------------------------------------------------------------------------

%% The rest is just plotting, what follow here is going to be very dependent on the particular situation
% and is just kept here as an example.

hold on ; plot(xgl/1000,ygl/1000,LineWidth=2,Color="r")
plot(xglc/1000,yglc/1000,"ok");


GLNormals=FindOrCreateFigure("GLNormals");

quiver(GLQ.xglc/1000,GLQ.yglc/1000,GLQ.nxGL,GLQ.nyGL);
hold on
plot(GLQ.xglc/1000,GLQ.yglc/1000,"-ok");
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF) ;
axis equal


%%
TN=FindOrCreateFigure("ThetaN");  clf(TN)  ; 
fprintf(' Plotting normal and tangential buttressing ratio \n ')


ax1=axes ;  % first axes
UaPlots(CtrlVar,MUA,F,"-speed-") ; %  plots (basal) speed
title(ax1,"")

% map=othercolor("YlOrRd9",1024) ;  colormap(ax1,map);
map=othercolor("BuPu4",1024) ;  colormap(ax1,map);
xlabel("xps (km)",Interpreter="latex")
ylabel("yps (km)",Interpreter="latex")

ax2=axes ; % second axes

% ACHTUNG:  Here is get rid of "outliers", this is possibly questionable, and should ony be done for plotting purposes and
% after having having checked that this is justified. 
I=isoutlier(GLQ.ThetaN) ;


x=GLQ.xglc(~I)/1000 ; y=GLQ.yglc(~I)/1000 ; Values=GLQ.ThetaN(~I) ;


% map=othercolor("YlOrRd9",Ncol) ;
Ncol=1024;
map=jet(Ncol);

% ACHTUNG: Now be carefull as here I truncate the values above and below given max and min values!!!
ValuesMin=-0.5 ;  ValuesMax=1.5 ;  Values(Values>ValuesMax)= ValuesMax ; Values(Values<ValuesMin)=ValuesMin ;
% ValuesMin=min(ValuesMax)  ; ValuesMax=max(ValuesMin) ; % better to do this first, and then possibly truncate afterwards

ind=round((Ncol-1)*(Values-ValuesMin)/(ValuesMax-ValuesMin))+1 ; % index into colormap
ValueOne=1;

indOneUpper=round((Ncol-1)*(ValueOne*1.1-ValuesMin)/(ValuesMax-ValuesMin))+1 ; % index into colormap
indOneLower=round((Ncol-1)*(ValueOne*0.9-ValuesMin)/(ValuesMax-ValuesMin))+1 ; % index into colormap

map(indOneLower:indOneUpper,:)=map(indOneLower:indOneUpper,:)*0 ;
% Note: For this colormapping to be correct, the limits of caxis must be manually set to [ValuesMin ValuesMax]
%
colormap(ax2,map);

hold on
% scatter(x,y,100*Values,map(ind,:),'filled');
bc=bubblechart(ax2,x,y,Values*0+1,map(ind,:)) ; %,'filled');
axis equal ;
hold on
plot(xglc/1000,yglc/1000,"-r");



xlabel("xps (km)",Interpreter="latex")
ylabel("yps (km)",Interpreter="latex")


linkaxes([ax1,ax2])
caxis(ax2,[min(Values),max(Values)])
hold off
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.07 .2 .05 .6]);
cb2 = colorbar(ax2,'Position',[.88 .2 .05 .6]);  caxis(ax2,[ValuesMin ValuesMax]) ; % 

title(cb2,"$\theta_N$",interpreter="latex")
title(cb1,"Speed (m/yr)",interpreter="latex")
%colormap(ax1,'hot')
%colormap(ax2,'cool')
% axis([-1700 -1100 -700 -200])
axis([-1700 -1400 -700 -200])
bubblesize([3 20])
TN.Position=[50 150 800 1100] ;
cb2.Position=[0.900 0.550 0.0300 0.3500];
cb1.Position=[0.900 0.1400 0.0300 0.3500];



% bl=bubblelegend("$\theta$",Interpreter="latex");

%%

Fig=gca ;

Region="Getz";
Region="Whole";

switch Region

    case "Whole"
        bubblesize([3 10])
        PlotLatLonGrid(1000,5/4,10/2); % All
        % exportgraphics(Fig,'C:\Users\lapnjc6\OneDrive - Northumbria University - Production Azure AD\Work\Manuscripts\2022 ThwaitesIceShelfButtressing\Figures\ButtressingRatioWhole.pdf')

    case "PIG"

        Fig.XLim=[-1680 -1550];  Fig.YLim=[-400 -200]; bl.Location="NorthWest";  % PIG
        % exportgraphics(Fig,'C:\Users\lapnjc6\OneDrive - Northumbria University - Production Azure AD\Work\Manuscripts\2022 ThwaitesIceShelfButtressing\Figures\ButtressingRatioPIG.pdf')

    case "Thwaites"

        Fig.XLim=[-1580 -1515];  Fig.YLim=[-480 -380]; bl.Location="NorthEast"; PlotLatLonGrid(1000,5/2^4,10/2^3); % Thwaites
        % exportgraphics(Fig,'C:\Users\lapnjc6\OneDrive - Northumbria University - Production Azure AD\Work\Manuscripts\2022 ThwaitesIceShelfButtressing\Figures\ButtressingRatioThwaites.pdf')


    case "Getz" % actually Pope, Smith and Kohler

        Fig.XLim=[-1600 -1480];  Fig.YLim=[-700 -520]; bl.Location="NorthWest";

    otherwise

        error("what case")

end

%
% exportgraphics(Fig,'C:\Users\lapnjc6\OneDrive - Northumbria University - Production Azure AD\Work\Manuscripts\2022 ThwaitesIceShelfButtressing\Figures\ButtressingRatioWhole.pdf')

return ; % but stuff below possibly also of interest...

% q= (A (rho g)^(n+1) (1-rho/rhow)^n' / (4^n' C') ))^(1/(1+m)) theta^(n/(m'+1)) h^((m'+n+3)/(m+1))  ( Schoof notation)
%
% where tau_bx=C' |\mathbf{v}|^{m'-1} u and |tau| = C' |v|^m' (Schoof)
% I use: u=C |tau|^{m-1} tau_bx
% and I have |v|= C |tau|^m
% therfore:  |v|= C (C' |v|^m')^m = C C'^m  |v|^(m' m)
%->  m'=1/m and C'=(1/C)^(1/m)
%
% % q= A (rho g)^(n+1) (1-rho/rhow)^n / (4^n C^(-1/m)) ) theta^(n/(1/m+1)) h^((1/m+n+3)/(1/m+1))  ( my notation)
%
% n=3 ; m=3;
% AGlen=1e-14*1e9*365.25*24*60*60; C0=7.624e6^(-m)*1000^m*365.25*24*60*60;
% rho=900 ; rhow=1000; g=9.81/1000;
lw=10;
A0=mean(AGlen) ; C0=mean(C);
qSchoof=rhogl.* (A0*(rhogl*g).^(n+1).*(1-rhogl/rhow).^n/(4^n*C0^(-1/m))).^(1/(1/m+1)).*ThetaN.^(n/(1/m+1)).*hgl.^((1/m+n+3)/(1/m+1)) ; % my notation
qSchoof=real(qSchoof);



%%



T=(qngl-qSchoof)./qngl;

figure
CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
figure ; set(gcf,'Position',figPosition)
PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ;
hold on;
color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,qngl*0+100,T,'LineWidth',lw) ; colorbar
set(gcf,'Name','(q-qSchoof)/q')   ;
colorbar ; xlabel('x (km)') ; ylabel('y (km)')
axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])

%%

fprintf(' Plotting ice flux normal and tangential to grounding line \n ')

CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
figure ; set(gcf,'Position',figPosition)
PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ;
hold on; htag=color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,qngl*0+100,qngl/1e9,'LineWidth',lw) ;
%title('Ice flux normal to grounding line (10^9 kg m^{-1} a^{-1})')   ;
set(gcf,'Name','Ice flux normal to grounding line (10^9 kg m^{-1} a^{-1})')   ;
colorbar ; xlabel('x (km)') ; ylabel('y (km)')
axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])



CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
figure ; set(gcf,'Position',figPosition)
PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ;
hold on; htag=color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,qtgl*0+100,qtgl/1e9,'LineWidth',lw) ;
title('Ice flux tangential to grounding line (10^9 kg m^{-1} a^{-1})')   ;
set(gcf,'Name','Ice flux tangential to grounding line (10^9 kg m^{-1} a^{-1})')   ;
colorbar ; xlabel('x (km)') ; ylabel('y (km)')
axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])


%%

fprintf(' Plotting normal and tangential buttressing ratio \n ')


CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
figure ; set(gcf,'Position',figPosition)
PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ; daspect([1 1 1]) ;
hold on; htag=color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,kappaN*0+100,kappaN,'LineWidth',lw) ;

set(gcf,'Name','\kappa_n')   ;
colorbar ; xlabel('x (km)') ; ylabel('y (km)')
axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])

CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
figure ; set(gcf,'Position',figPosition)
PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ; daspect([1 1 1]) ;
hold on; color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,kappaT*0+100,kappaT,'LineWidth',lw) ;
%title('\Theta_t')   ;
set(gcf,'Name','\kappa_t')   ;
colorbar ; xlabel('x (km)') ; ylabel('y (km)')
axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])


fprintf(' Plotting deviatoric stresses along grounding line \n ')
CtrlVar.PlotNodes=1; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
figure ; set(gcf,'Position',figPosition)
%PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ; daspect([1 1 1]) ;

% grounding line profile
%plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'g','LineWidth',lw)
plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'g','LineWidth',2)
%plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'g','LineWidth',2)

hold on
scale=0.05 ;

ind=1:50:numel(xglc);

% label
xls=(fxmin+0.1*(fxmax-fxmin))*CtrlVar.PlotXYscale; yls=(fymax-0.1*(fymax-fymin))*CtrlVar.PlotXYscale;
xlst=(fxmin+0.2*(fxmax-fxmin))*CtrlVar.PlotXYscale; ylst=(fymax-0.1*(fymax-fymin))*CtrlVar.PlotXYscale;
tx=100 ; ty=-100;
text(xlst/CtrlVar.PlotXYscale,ylst/CtrlVar.PlotXYscale,'(100 kPa)','HorizontalAlignment','center')

xlu=(fxmin+0.07*(fxmax-fxmin))*CtrlVar.PlotXYscale; ylu=(fymax-0.2*(fymax-fymin))*CtrlVar.PlotXYscale;
xlut=(fxmin+0.2*(fxmax-fxmin))*CtrlVar.PlotXYscale; ylut=(fymax-0.2*(fymax-fymin))*CtrlVar.PlotXYscale;
us=500 ; vs=0;
text(xlut/CtrlVar.PlotXYscale,ylut/CtrlVar.PlotXYscale,'(500 m/a)','HorizontalAlignment','center')


plot_strains([xls ;xglc(ind)]/CtrlVar.PlotXYscale,[yls;yglc(ind)]/CtrlVar.PlotXYscale,[tx;txxgl(ind)],[0;txygl(ind)],[ty;tyygl(ind)],scale,1);
quiver([xlu;xglc(ind)]/CtrlVar.PlotXYscale,[ylu;yglc(ind)]/CtrlVar.PlotXYscale,[us;ugl(ind)],[vs;vgl(ind)],'color','k')
axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])



%title('Deviatoric stresses') ;
xlabel('x (km)') ; ylabel('y (km)')
set(gcf,'Name','Deviatoric stresses')   ;
slope0a=973.6686e3;  slope0b=1265.712e3;
y=coordinates(:,2);
hold on ; plot([slope0a/CtrlVar.PlotXYscale,slope0a/CtrlVar.PlotXYscale],[min(y) max(y)]/CtrlVar.PlotXYscale,'m')
hold on ; plot([slope0b/CtrlVar.PlotXYscale,slope0b/CtrlVar.PlotXYscale],[min(y) max(y)]/CtrlVar.PlotXYscale,'m')

%%

end

%%
