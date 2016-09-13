function [UserVar,ub1,vb1,h1,luv1,lh1,RunInfo]=...
    SSTREAM_TransientImplicit(UserVar,CtrlVar,MUA,BCs,dt,h0,S,B,ub0,vb0,ub1,vb1,h1,as0,ab0,as1,ab1,dudt,dvdt,luv,lh,...
    AGlen,C,n,m,alpha,rho,rhow,g)



nOut=nargout;
if nOut~=7
    error('Ua:SSTREAM_TransientImplicit','Need 7 output arguments')
end


if CtrlVar.InfoLevelNonLinIt>=10  ; fprintf(CtrlVar.fidlog,' \n SSTREAM(uvh): Transient implicit with respect to u, v, and h  \n ') ; end

%%
% Fully implicit Newton-Raphson with regard to both u, v and h
% advances the solution by dt
%
% h0, u0, and v0 are values at the start of the time step,
% on input h1,u1,v1 are estimates for h, u, and v at the end of the time step
% on exit  h1, u1, and v1 are calculated values for u,v and h at the end of the time step



%%



% I need to solve
%
% [Kxu Kxv Kxh Luv'  0  ] [du]        =  [ -Ru ] - Luv' luv
% [Kyu Kyv Kyh          ] [dv]           [ -Rv ]
% [Khu Khv Khh  0   Lh' ] [dh]           [ -Rh- Lh' lh ]
% [  Luv        0    0  ] [duv]          [ cuv-Luv [u ;v ]
% [ 0    Lh  0  0    0  ] [dlh]          [ ch-Lh h]
%
% All matrices are Nnodes x Nnodes, apart from:
% Luv is #uv constraints x 2 Nnodes, i.e. Luv [u;v]= cuv
% Lh  is # h contraints x Nnodes, i.e.    Lh h= ch
%  or
%
% [K L'] [ duvh ]      =  [ -R- L' l ]
% [L 0 ] [  dl  ]         [ cuvh-L [u;v;h]  ]
%
% where
%
% K= [Kxu Kxv Kxh]
%    [Kyu Kyv Kyh]
%    [Khu Khv Khh]
%
% and
% L=[Luv 0]
%   [0  Lh]
% and uvh=[u;v;h], duvh=[du;dv; dh]  and l=[luv ; lh]
% where L [u;v;h]=cuvh
% 
 
            
if any(h0<0) ; warning('MATLAB:SSTREAM_TransientImplicit',' thickness negative ') ; end


tStart=tic;



ub=ub1 ; vb=vb1 ; h=h1 ; % starting values
dub=ub-ub0; dvb=vb-vb0 ; dh=h-h0;
%du=zeros(Nnodes,1) ; dv=zeros(Nnodes,1) ; dh=zeros(Nnodes,1) ;


%% assemble global Lagrange constraint matrix
MLC=BCs2MLC(MUA,BCs);
Luv=MLC.ubvbL;  
cuv=MLC.ubvbRhs;
Lh=MLC.hL;
ch=MLC.hRhs;

if numel(luv)~=numel(cuv) ; luv=zeros(numel(cuv),1) ; end
if numel(lh)~=numel(ch) ; lh=zeros(numel(ch),1) ; end
nluv=numel(luv) ; nlh=numel(lh); 

[L,cuvh,l]=AssembleLuvh(Luv,Lh,cuv,ch,luv,lh,MUA.Nnodes);
dl=l*0;


%% calculate basis function derivatives
% this is a very memory expensive way of doing this, but it saves
% redoing this calculation within every non-lin iteration

% Deriv : Nele x dof x nod
%  detJ : Nele
%dof=2; [Nele,nod]=size(connectivity);
%MeshProp.Deriv=zeros(Nele,dof,nod,nip);
%MeshProp.DetJ=zeros(Nele,nip);



iteration=0 ; Stagnated=0;
r=1e10; diffVector=zeros(CtrlVar.NRitmax,1); diffDu=1e10 ; diffDh=1e10; diffDlambda=1e10;

while (r>1e-15 && ((r> CtrlVar.NLtol || diffDu > CtrlVar.du || diffDh> CtrlVar.dh  || diffDlambda > CtrlVar.dl) &&  iteration <= CtrlVar.NRitmax && ~Stagnated) ) || iteration < CtrlVar.NRitmin
    iteration=iteration+1;
    
    % If I want to implement an in complete NR method, then presumably all that is
    % needed is to call uvhAssembly once in while with just the first output
    % argument, ie R, and only update K occationally
    [UserVar,R,K,~,FI]=uvhAssembly(UserVar,CtrlVar,MUA,ub,vb,h,S,B,ub0,vb0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g);
    
    if iteration==1 
        F0=FI ; % F0 is used as a normalisation factor when calculating the residual, do not change this normalisation factor in the course of the iteration
    end  
    
    
    %% solve the (asymmetrical) linear system
    if ~isempty(L)
        frhs=-R-L'*l;
        grhs=cuvh-L*[ub;vb;h];
    else
        frhs=-R;
        grhs=[];
    end
    
   [r0,rl0,ruv0,rh0]=ResidualCostFunction(frhs,grhs,F0,MUA.Nnodes); 

    [duvh,dl]=solveKApe(K,L,frhs,grhs,[dub;dvb;dh],dl,CtrlVar);
    
    
    if any(isnan(duvh)) ; save TestSave  ; fprintf(CtrlVar.fidlog,'error: NaN in solution of implicit system \n') ; error(' NaN in solution of implicit uvh system ' ) ; end
    %%
    dub=duvh(1:MUA.Nnodes) ;  dvb=duvh(MUA.Nnodes+1:2*MUA.Nnodes); dh=duvh(2*MUA.Nnodes+1:end);
    
    
    %% calculate  residuals at full Newton step
    gamma=1;
    [UserVar,r1,ruv1,rh1,rl1]=CalcCostFunctionNRuvh(UserVar,CtrlVar,MUA,gamma,dub,dvb,dh,ub,vb,h,S,B,ub0,vb0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g,F0,L,l,dl,cuvh);
                          
    
    %% either accept full Newton step or do a line search
    
    [UserVar,r,ruv,rh,rl,gamma,infovector,iarm,BacktrackInfo]=FindBestGamma2DuvhBacktrack...
        (UserVar,CtrlVar,MUA,F0,r0,r1,ruv1,rh1,rl1,ub,vb,h,dub,dvb,dh,S,B,ub0,vb0,h0,L,l,dl,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g,cuvh);
    
    
    if BacktrackInfo.converged==0
        Stagnated=1;
    end
    
    %% If desired, plot residual along search direction
    if CtrlVar.InfoLevelNonLinIt>=10 && CtrlVar.doplots==1
        nnn=12;
        gammaTestVector=zeros(nnn,1) ; rTestvector=zeros(nnn,1);
        Up=2.2;
        if gamma>0.7*Up ; Up=2*gamma; end
        parfor I=1:nnn
            gammaTest=Up*(I-1)/(nnn-1)+gamma/50;
            [~,rTest,ruv1Test,rh1Test,rl1Test]=CalcCostFunctionNRuvh(UserVar,CtrlVar,MUA,gammaTest,dub,dvb,dh,ub,vb,h,S,B,ub0,vb0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g,F0,L,l,dl,cuvh);
            gammaTestVector(I)=gammaTest ; rTestvector(I)=rTest;
        end
        gammaTestVector=[gammaTestVector(:);infovector(:,1)];
        rTestvector=[rTestvector(:);infovector(:,2)];
        [gammaTestVector,ind]=unique(gammaTestVector) ; rTestvector=rTestvector(ind) ;
        [gammaTestVector,ind]=sort(gammaTestVector) ; rTestvector=rTestvector(ind) ;
        
        
        slope=-2*r0;
        figure ; plot(gammaTestVector,rTestvector,'o-r') ; hold on ;
        plot(gamma,r,'Marker','h','MarkerEdgeColor','k','MarkerFaceColor','g')
        plot([gammaTestVector(1) gammaTestVector(2)],[rTestvector(1) rTestvector(1)+(gammaTestVector(2)-gammaTestVector(1))*slope],'g')
        
        title(sprintf('uvh iteration %-i,  iarm=%-i ',iteration,iarm)) ; xlabel(' \gamma ') ; ylabel('Residual')
        
        
        
        
        
        hold off
        %input('press return to continue')
    end
    
    
    %% calculate statistics on change in speed, thickness and Lagrange parameters
    D=mean(sqrt(ub.*ub+vb.*vb))+CtrlVar.SpeedZero;
    diffDu=gamma*full(max(abs(dub))+max(abs(dvb)))/D;        % sum of max change in du and dv normalized by mean speed
    diffDh=gamma*full(max(abs(dh))/(mean(abs(h)))+0.01);            % max change in thickness divided by mean thickness
    diffDlambda=gamma*max(abs(dl))/mean(abs(l));
    diffVector(iteration)=r0;   % override last value, because it was just an (very accurate) estimate
    diffVector(iteration+1)=r;
    
    
    if isempty(diffDlambda)
        diffDlambda=0;
    end
    
    %% update variables
    ub=ub+gamma*dub ; vb=vb+gamma*dvb; h=h+gamma*dh;  l=l+gamma*dl;
    
    luv=l(1:nluv) ; lh=l(nluv+1:end);
    
    if~isempty(Lh)
        BCsNormh=norm(ch-Lh*h);
    else
        BCsNormh=0;
    end
    
    if~isempty(Luv)
        BCsNormuv=norm(cuv-Luv*[ub;vb]);
    else
        BCsNormuv=0;
    end
    
    
    temp=CtrlVar.ResetThicknessToMinThickness;
    if ~CtrlVar.ResetThicknessInNonLinLoop
        CtrlVar.ResetThicknessToMinThickness=0;
    end
    [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
    CtrlVar.ResetThicknessToMinThickness=temp;
    
    if CtrlVar.MassBalanceGeometryFeedback==1
        GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
        rdamp=CtrlVar.MassBalanceGeometryFeedbackDamping;
        if rdamp~=0
            as1Old=as1 ; ab1Old=ab1;
        end
        [as1,ab1]=GetMassBalance(CtrlVar.UserVar,CtrlVar,MUA,CtrlVar.time+dt,s,b,h,S,B,rho,rhow,GF);
        
        if rdamp~=0
            % If Hessian inaccurate, or too non-linear, then
            % dampen these changes might be a good idea.
            as1=(1-rdamp)*as1+rdamp*as1Old;
            ab1=(1-rdamp)*ab1+rdamp*ab1Old;
        end
    end

    
    if CtrlVar.InfoLevelNonLinIt>=100  && CtrlVar.doplots==1
        PlotForceResidualVectors('uvh',R,L,l,MUA.coordinates,CtrlVar) ; axis equal tight
    end
    
    if CtrlVar.InfoLevelNonLinIt>=1
        fprintf(CtrlVar.fidlog,'NR-STREAM(uvh):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , ruv=%-14.7g , rh=%-14.7g , du=%-14.7g , dh=%-14.7g , dl=%-14.7g , BCsNormuv=%-g , BCsNormh=%-g  \n ',...
            iteration,iarm,gamma,r/r0,r0,r,ruv,rh,diffDu,diffDh,diffDlambda,BCsNormuv,BCsNormh);
    end
    
end

%% return calculated values at the end of the time step
ub1=ub ; vb1=vb ; h1=h; luv1=luv ; lh1=lh;

tEnd=toc(tStart);

%% print/plot some info

if CtrlVar.InfoLevelNonLinIt>=10 && iteration >= 2 && CtrlVar.doplots==1
    
    N=max([1,iteration-5]);
    
    [~,~,a1]=detrend_xt(log10(diffVector(N:iteration)),N:iteration);
    fprintf(CtrlVar.fidlog,' slope NR : %14.7g \n',a1);
    figure; semilogy(0:iteration,diffVector(1:iteration+1),'x-r') ; title('NR uvh implicit') ; xlabel('Iteration') ; ylabel('Residual')
end

if ~isempty(L)
    BCerror=norm(L*[ub;vb;h]-cuvh);
    if BCerror>0
        fprintf(CtrlVar.fidlog,'Error in satisfying Dirichlet BC %14.7g  \n ',norm(L*[ub;vb;h]-cuvh));
    end
end


if r>CtrlVar.NLtol
    fprintf(CtrlVar.fidlog,' SSTREAM(uvh)  did not converge to given tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',CtrlVar.NLtol,r,iteration,tEnd);
    fprintf(CtrlVar.InfoFile,' SSTREAM(uvh)  did not converge to given tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',CtrlVar.NLtol,r,iteration,tEnd);
    
    RunInfo.converged=0;
else
    RunInfo.converged=1;
    if CtrlVar.InfoLevelNonLinIt>=1
        fprintf(CtrlVar.fidlog,' SSTREAM(uvh) (time|dt)=(%g|%g): Converged to given tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',...
            CtrlVar.time,dt,CtrlVar.NLtol,r,iteration,tEnd) ;
        
        fprintf(CtrlVar.InfoFile,' SSTREAM(uvh/%s) \t time=%15.5f \t dt=%-g \t r=%-g \t #it=% i \t CPUsec=%-g \n',...
            CtrlVar.uvhTimeSteppingMethod,CtrlVar.time,dt,r,iteration,tEnd) ;
    end
end


if iteration > CtrlVar.NRitmax
    fprintf(CtrlVar.fidlog,'Warning: maximum number of NRuvh iterations %-i reached \n',CtrlVar.NRitmax);
    warning('SSTREAM2dNR:MaxIterationReached','SSTREAM2NR exits because maximum number of iterations %-i reached \n',CtrlVar.NRitmax)
end

RunInfo.Iterations=iteration;   RunInfo.residual=r;



%     if any(isnan(u1)) ; save TestSave  ;  error(' NaN in u1 ' ) ; end
%     if any(isnan(v1)) ; save TestSave  ;  error(' NaN in v1 ' ) ; end
%     if any(isnan(h1)) ; save TestSave  ;  error(' NaN in h1 ' ) ; end
%     if any(isnan(lambdauv)) ; save TestSave  ;  error(' NaN in lambdauv ' ) ; end
%     if any(isnan(lambdah)) ; save TestSave  ;  error(' NaN in lambdah ' ) ; end


end


