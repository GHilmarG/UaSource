function ab=JenkinsIceShelfMeltParameteristation(CtrlVar,MUA,s,b,S,B,GF,Tw,Sw,tcDe)


%fprintf(' Jenkins \n ')
[dfdx,dfdy]=calcFEderivativesMUA(b,MUA,CtrlVar);
bGrad=sqrt(dfdx.*dfdx+dfdy.*dfdy);  % defined at integration points
[M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
bGrad=ProjectFintOntoNodes(MUA,bGrad);
%bGrad=M*bGrad;
%         for I=1:10
%             bGrad=Nodes2EleMean(MUA.connectivity,bGrad);
%             bGrad=M*bGrad;
%         end

minSlope=0.001;
bGrad(bGrad<minSlope)=minSlope;


% find for each node on an ice shelf the min distance to a grounded point that has
% a greater (more negative) draft than the point itself
jf=find(GF.node<0.5) ;  % These nodes are afloat
jg=FindNearestGroundedPoint(CtrlVar,MUA,s,b,S,B,GF,jf) ;

glDe=b;
glDe(jf)=b(jg);

%Tw=CtrlVar.Tw+CtrlVar.DeltaTw*sin(2*pi*time/CtrlVar.PeriodTw);
%Tw=CtrlVar.Tw;
%Sw=CtrlVar.Sw;
%tcDe=CtrlVar.tcDe;
%fprintf('Jenkins melt rate parameters: Tw=%f \t Sw=%f \t glDe=%f \n',Tw,Sw,glDe)

ab(jf)=basal_melt(b(jf),bGrad(jf),glDe(jf),Tw,Sw,tcDe);


ab=-ab; % Jenkins defines pos as melting

% keep within limits
[ab,iU,iL] = kk_proj(ab,500,-500);


end