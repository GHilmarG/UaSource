
function [exx,eyy,exy,e]=CalcNodalStrainRates(CtrlVar,MUA,u,v)


%%
%  
% Calculates horizontal strain rates given horizontal velocities. 
%
%   [exx,eyy,exy,e]=CalcNodalStrainRates(CtrlVar,MUA,u,v)
%
% Returns nodal values.
%
% Note: The strain rates are calculated at integration points and then projected
% onto nodes. 
%
% The projection does not conserve positivity and positve integration
% values can become negative at nodes. The effectiv strain rate , e, is for
% this reason calculated directly from nodal values, ensuring that e is
% always positive.
%
% Example: 
%
% Read data, calculate nodal strain rates, and then plot over FE mesh at roughly equal spaced grid.
%
%   load ('GaussPeak_Example_Restartfile.mat','MUA','CtrlVarInRestartFile','F','GF','BCs');  % load data
%   CtrlVar=CtrlVarInRestartFile; x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2);
%   [exx,eyy,exy,e]=CalcNodalStrainRates(CtrlVar,MUA,F.ub,F.vb);                             % calculate strain rates
%   [X,Y]=ndgrid(linspace(min(x),max(x),10),linspace(min(y),max(y),10));
%   I=nearestNeighbor(MUA.TR,[X(:) Y(:)]);  % find nodes within computational grid closest to the regularly scape X and Y grid points.
%   figure
%   CtrlVar.PlotNodes=0; PlotMuaMesh(CtrlVar,MUA,[],'color','k') ;                           % Plot FE mesh
%   hold on
%   scale=1e5; LineWidth=2 ;
%   PlotTensor(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,exx(I),exy(I),eyy(I),scale,LineWidth);  % plot strain rates
%   axis equal
%
%
% See also CalcNodalStrainRatesAndStresses


unod=reshape(u(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(v(MUA.connectivity,1),MUA.Nele,MUA.nod);
exx=zeros(MUA.Nele,MUA.nip);
eyy=zeros(MUA.Nele,MUA.nip);
exy=zeros(MUA.Nele,MUA.nip);

for Iint=1:MUA.nip
    
    Deriv=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    
    for I=1:MUA.nod
        exx(:,Iint)=exx(:,Iint)+Deriv(:,1,I).*unod(:,I);
        eyy(:,Iint)=eyy(:,Iint)+Deriv(:,2,I).*vnod(:,I);
        exy(:,Iint)=exy(:,Iint)+0.5*(Deriv(:,1,I).*vnod(:,I) + Deriv(:,2,I).*unod(:,I));
    end
end


[exx,eyy,exy]=ProjectFintOntoNodes(MUA,exx,eyy,exy);
 
e=real(sqrt(exx.^2+eyy.^2+exx.*eyy+exy.^2));

end
