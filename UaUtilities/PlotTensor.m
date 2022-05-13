function PlotTensor(x,y,txx,txy,tyy,scale,LineWidth)

%%
% Plots a 2x2 symmetrical tensor quantity
%
%   PlotTensor(x,y,txx,txy,tyy,scale,LineWidth)
%
%   T=[txx  txy]
%     [txy  tyy]
%
% Compression is plotted in red, extension in blue.
%
% *Example:*
%
%   load ('GaussPeak_Example_Restartfile.mat','MUA','CtrlVarInRestartFile','F','GF','BCs');  % load data
%   CtrlVar=CtrlVarInRestartFile; x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2);
%   [exx,eyy,exy,e]=CalcNodalStrainRates(CtrlVar,MUA,F.ub,F.vb);                             % calculate strain rates
%   [X,Y]=ndgrid(linspace(min(x),max(x),20),linspace(min(y),max(y),20));
%   I=nearestNeighbor(MUA.TR,[X(:) Y(:)]);  % find nodes within computational grid closest to the regularly scape X and Y grid points.
%   FigTensor=figure;
%   CtrlVar.PlotNodes=0; PlotMuaMesh(CtrlVar,MUA,[],'color','k') ;                           % Plot FE mesh
%   hold on
%   scale=1e4; LineWidth=2 ;
%   PlotTensor(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,exx(I),exy(I),eyy(I),scale,LineWidth);  % plot strain rates
%   axis equal tight
%   FigTensor.Children.Title.String='Strain rates'
%
%%


% figure ; PlotTensor(F.x/1000,F.y/1000,F.exx,F.eyy,F.exy,1) ; axis equal
if nargin<7 
    LineWidth=1; 
end

headscale=0.3;
sharp=0.3;



for k=1:length(x)
    
    if ~isnan(txx(k))
        D=[txx(k) txy(k) ; txy(k) tyy(k)];
        [pAxis,pStrains]=eig(D); l1=pStrains(1,1) ; l2=pStrains(2,2);
        
        % pStrains : Principla strains
        
        p1x=l1*pAxis(1,1) ; p1y=l1*pAxis(2,1) ;
        p2x=l2*pAxis(1,2) ; p2y=l2*pAxis(2,2) ;
        
        
        if l1 < 0
            head=-1;
            col=[1 0 0 ];
        else
            head=1;
            col=[0 0 1];
        end
        
        
        ghg_arrow(x(k),y(k),p1x,p1y,scale,headscale,sharp,head,col,LineWidth);
        ghg_arrow(x(k),y(k),-p1x,-p1y,scale,headscale,sharp,head,col,LineWidth);
        
        if l2 < 0
            head=-1;
            col=[1 0 0 ];
        else
            head=1;
            col=[0 0 1];
        end
        ghg_arrow(x(k),y(k),p2x,p2y,scale,headscale,sharp,head,col,LineWidth);
        ghg_arrow(x(k),y(k),-p2x,-p2y,scale,headscale,sharp,head,col,LineWidth);
    end
end


return

end
