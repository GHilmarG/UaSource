function [UserVar,dJdC,dJdAGlen,ub,vb,ud,vd,uAdjoint,vAdjoint,dIdCreg,dIdAGlenreg,dIdCdata,dIdAGlendata,dIdCbarrier,dIdAGlenbarrier,lambdaAdjoint]=...
    AdjointGradientNR2d(...
    UserVar,CtrlVar,MUA,BCs,BCsAdjoint,s,b,h,S,B,ub,vb,ud,vd,uvAdjoint,AGlen,C,n,m,alpha,rho,rhow,g,GF,Priors,Meas)

nargoutchk(16,16)
narginchk(26,26)

if CtrlVar.AGlenisElementBased
    dIdAGlendata=zeros(MUA.Nele,1) ;
else
    dIdAGlendata=zeros(MUA.Nnodes,1);
end

if CtrlVar.CisElementBased
    dIdCdata=zeros(MUA.Nele,1);
else
    dIdCdata=zeros(MUA.Nnodes,1);
end




%% Step 1: solve linearized forward problem
%[ub,vb,ud,vd,l,Kuv,Ruv,RunInfo]= uv(CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);

[UserVar,ub,vb,ud,vd,uvAdjoint,Kuv,Ruv,RunInfo,ubvbL]=uv(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,uvAdjoint,AGlen,C,n,m,alpha,rho,rhow,g,GF);

%% Step 2:  Solve adjoint equation, i.e.   K l=-r
% fprintf(' Solve ajoint problem \n ')
% I need to impose boundary conditions on lx and ly
% if the problem is (fully) adjoint I have exactly the same BC
% I need to solve
%
% [Kxu Kxv Luv'] [lx]        =  [ u-uMeas ]
% [Kyu Kyv     ] [ly]        =  [ v-vMeas ]
% [  Luv      0] [lambdauv]     [ Luvrhs  ]
% All matrices are Nnodes x Nnodes, apart from:
% Luv is #uv constraints x 2 Nnodes


% forming right-hand side of the adjoint equations

[J,Idata,IRegC,IRegAGlen,dIduv,IBarrierC,IBarrierAGlen]=MisfitFunction(UserVar,CtrlVar,MUA,ub,vb,ud,vd,AGlen,C,Priors,Meas);


rhs=dIduv(:);


MLC_Adjoint=BCs2MLC(MUA,BCsAdjoint);
LAdjoint=MLC_Adjoint.ubvbL;
LAdjointrhs=MLC_Adjoint.ubvbRhs;
lambdaAdjoint=zeros(numel(LAdjointrhs),1) ;

[uvAdjoint,lambdaAdjoint]=solveKApeSymmetric(Kuv,LAdjoint,rhs,LAdjointrhs,[],lambdaAdjoint,CtrlVar);


if ~isreal(uvAdjoint) 
    save TestSave ; error('When solving adjoint equation Lagrange parmeters complex ')
end

uAdjoint=real(uvAdjoint(1:MUA.Nnodes)) ; vAdjoint=real(uvAdjoint(MUA.Nnodes+1:2*MUA.Nnodes));

if CtrlVar.InfoLevelAdjoint>=1000 && CtrlVar.doplots
    
    GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
    tri=MUA.connectivity;
    
    figure
    hold off
    subplot(2,2,1)
    [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,dIduv(1:length(ub)),CtrlVar);  title('dIdu')
    hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        
    subplot(2,2,2)
    [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,dIduv(1+length(ub):end),CtrlVar);  title('dIdv')
    hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
    
    subplot(2,2,3)
    [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,uAdjoint,CtrlVar);  title('lx')
    hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
    
    subplot(2,2,4)
    [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,vAdjoint,CtrlVar);  title('ly')
    hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
end


switch upper(CtrlVar.AdjointGrad)
    
    case 'C'
        
        switch lower(CtrlVar.AdjointGradientEvaluation)
            
            case 'discrete' % Direct gradient evaluated from nodal points.
                
                if CtrlVar.CisElementBased
                    M= Ele2Nodes(MUA.connectivity,MUA.Nnodes);
                    Cnode=M*C;
                else
                    Cnode=C;
                end
                
                dIdCdata = -(1/m)*GF.node.*(Cnode+CtrlVar.CAdjointZero).^(-1/m-1).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1/m-1).*(u.*uAdjoint+v.*vAdjoint);
                
                if CtrlVar.CisElementBased
                    dIdCdata=Nodes2EleMean(MUA.connectivity,dIdCdata);
                end
                
            case 'integral'
                
                if CtrlVar.CisElementBased
                    
                    dIdCdata=dIdCqEleSteps(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF);
                    
                else
                    dIdCdata=dIdCq(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF);
                end

                
            otherwise
                error(' what case ? ' )
        end
        
    case 'A'
        switch lower(CtrlVar.AdjointGradientEvaluation)
            
            case 'discrete' % Direct gradient evaluated from nodal points.
                
                
                fprintf(' CtrlVar.AdjointGradientEvaluation=''uvdiscrete'' not possible in a combination with AGlen inverstion\n')
                error('AdjointGradientNR2d:DiscreteAdjointAGlen','Discrete case not implemented. Used integral evaluation instead.')
            case 'integral'
                if CtrlVar.AGlenisElementBased
                    
                    dIdAGlendata=dIdAEleSteps(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF);
                                        
                else
                    
                    dIdAGlendata=dIdAq(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF);
                    
                end
            otherwise
                error(' what case ? ' )
        end
        
end



dIdCreg=Calc_dIregdC(CtrlVar,MUA,Priors.CovC,C,Priors.C);

%save TestSave ; error('fdsa')

dIdAGlenreg=Calc_dIregdAGlen(CtrlVar,MUA,Priors.CovAGlen,AGlen,Priors.AGlen);
dIdCbarrier=Calc_dIdCbarrier(CtrlVar,MUA,C);
dIdAGlenbarrier=Calc_dIdAGlenbarrier(CtrlVar,MUA,AGlen);


dIdCdata=CtrlVar.MisfitMultiplier*dIdCdata/CtrlVar.AdjointfScale;
dIdCreg=dIdCreg/CtrlVar.AdjointfScale;
dIdCbarrier=dIdCbarrier/CtrlVar.AdjointfScale;

dIdAGlendata=CtrlVar.MisfitMultiplier*dIdAGlendata/CtrlVar.AdjointfScale;
dIdAGlenreg=dIdAGlenreg/CtrlVar.AdjointfScale;
dIdAGlenbarrier=dIdAGlenbarrier/CtrlVar.AdjointfScale;

dJdC=dIdCdata+dIdCreg+dIdCbarrier;
dJdAGlen=dIdAGlendata+dIdAGlenreg+dIdAGlenbarrier;

if any(isnan(dJdC)) ; save TestSave ; error('NaN in dJdC') ; end
if any(isnan(dJdAGlen)) ; save TestSave ; error('NaN in dJdAGlen') ; end
if ~isreal(dJdC) ; save TestSave ; error('dJdC complex') ;end
if ~isreal(dJdAGlen) ; save TestSave ; error('dJdAGlen complex') ;end

%dJdC=real(dJdC);
%dJdAGlen=real(dJdAGlen);

if CtrlVar.InfoLevelAdjoint>=1000 && CtrlVar.doplots
    if ~isempty(strfind(CtrlVar.AdjointGrad,'C'))
        figure
        hold off
        subplot(2,2,1)
        %PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dIdCdata,CtrlVar) ; 
        
        PlotMeshScalarVariable(CtrlVar,MUA,dIdCdata);
        title('dIdCdata')
        hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        
        subplot(2,2,2)
        %PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dIdCreg,CtrlVar) ; 
        PlotMeshScalarVariable(CtrlVar,MUA,dIdCreg);
        title('dIregdC')
        hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        
        subplot(2,2,3)
        %PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dIdCbarrier,CtrlVar) ; 
        PlotMeshScalarVariable(CtrlVar,MUA,dIdCbarrier);
        title('dIdCbarrier')
        hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        
        subplot(2,2,4)
        %PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dJdC,CtrlVar) ;
        PlotMeshScalarVariable(CtrlVar,MUA,dJdC);
        title('dJdC')
        hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
    end
end

return

%% Testing gradient using brute force method

% calculates d kv/ d b using the complex number method
% using the fact that df/dx=Im(f(x+i dx))/dx
% dJ/dC=\p J/\p C + l^T [ \p f/\p C - \p(K u]/\p C]
%
% where K l =\p J /\p u
%
% Here K and u are `extended' states, i.e. including BCs, and K z = f
%
% So if J only depends on C implicity (i.e. \p J/\p C=0) then
%
% dJ/dC=\p ( l^t K z) /\p C,  where only the derivative of K is considerd
%

% 	rhs=Cd*[u-uMeas ; v-vMeas]; rhs=FErhs(rhs(1:Nnodes),rhs(Nnodes+1:end),coordinates,connectivity,nip);
% 	Luvrhs=Luvrhs*0;  % derivatives with respect to C are zero
%
% 	[l,lambdaAdjoint]=solveKApeSymmetric(K,Luv,rhs,Luvrhs,l,lambdaAdjoint,iteration,InfoLevelUzawa,SolMethod);
% 	%l=[l;lambdaAdjoint]; z=[u;v;lambdauv];

z=[u;v];
deltaC=1e-15;  Cpert=C;
%z=z*0+1 ; l=l*0+1;
%Nnorm=BasisFunctionNorms(connectivity,coordinates,nip);

dIdC=zeros(MUA.Nele,1);
for I=1:MUA.Nele
    %dIdC=zeros(Nnodes,1);
    Cpert(connectivity(I,:))=C(connectivity(I,:))+sqrt(-1)*deltaC; %/Nnorm(I);
    %Cpert(I)=C(I)+deltaC;
    [K,~,~,~]=KRTF(s,S,B,h,VUA.u,VUA.v,AGlen,n,Cpert,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar);
    Cpert(connectivity(I,:))=C(connectivity(I,:));
    dKdC=imag(K)/deltaC;
    %dKdC=(K-K0)/deltaC;
    dIdC(I)=l'*dKdC*z;
    
    %dIdC(I)=FE_inner_product(l,dKdC*z,coordinates,connectivity,nip);
    %figure(I+10000) ; trisurf(TRIxy,x,y,dIdC) ;  title('dIdC test ') ; xlabel('x') ; ylabel('y')
end

[xEle,yEle,DTele]=ElementCoordinates(connectivity,coordinates);

%figure(900) ; trisurf(DTele.Triangulation,xEle,yEle,dIdC) ;  title('dIdC ele ') ; xlabel('x') ; ylabel('y')

[dIdC]=dIdCqEleSteps(VUA,lx,ly,S,B,h,connectivity,coordinates,nip,C,m,rho,rhow,GF,CtrlVar)


%DkvDb=imag(kv)/db;
%DrhDb=imag(rh)/db;

end

%%




