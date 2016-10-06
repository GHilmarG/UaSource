function [R,K,F,T]=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt)

if nargout==1 ; OnlyR=1 ; else OnlyR=0 ; end

if CtrlVar.InfoLevelCPU>=10 ;   tCPU=tic; end

% NR:  h^{i+1} = \Delta h + h_1  and h^i=h_0
%      I solve for \Delta h.
%      If the thing converges then \Delta h -> 0 as i->\infty

h0=s0-b0;
h1=s1-b1;


ndim=2; neq=MUA.Nnodes;

theta=CtrlVar.theta;

h0nod=reshape(h0(MUA.connectivity,1),MUA.Nele,MUA.nod);
h1nod=reshape(h1(MUA.connectivity,1),MUA.Nele,MUA.nod);
s0nod=reshape(s0(MUA.connectivity,1),MUA.Nele,MUA.nod);
s1nod=reshape(s1(MUA.connectivity,1),MUA.Nele,MUA.nod);
a0nod=reshape(a0(MUA.connectivity,1),MUA.Nele,MUA.nod);
a1nod=reshape(a1(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~CtrlVar.AGlenisElementBased
    AGlen=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
    n=reshape(n(MUA.connectivity,1),MUA.Nele,MUA.nod);
end


% can I not just use h0 for h0MUA.nod, ie reuse variable?

[points,weights]=sample('triangle',MUA.nip,ndim);

K=sparse(neq,neq);

if ~OnlyR
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
else
    d1d1=[];
end

t1=zeros(MUA.Nele,MUA.nod); f1=zeros(MUA.Nele,MUA.nod);


% vector over all elements for each integartion point

if CtrlVar.ParallelAssembly
    
    parfor Iint=1:MUA.nip
        [tt,ff,dddd]=SSHEETintAssembly(Iint,CtrlVar,MUA,AGlen,n,rhonod,g,s0nod,h0nod,s1nod,h1nod,a0nod,a1nod,dt,OnlyR);
        t1=t1+tt;
        f1=f1+ff;
        d1d1=d1d1+dddd;
    end
    
else
    
    
    
    for Iint=1:MUA.nip
        
        %% [--
        fun=shape_fun(Iint,ndim,MUA.nod,points) ; % MUA.nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
        
        %     % Deriv : MUA.Nele x dof x MUA.nod
        %  detJ : MUA.Nele
        
        % values at integration point
        
        
        h0int=h0nod*fun;
        h1int=h1nod*fun;
        a0int=a0nod*fun;
        a1int=a1nod*fun;
        rhoint=rhonod*fun;
        
        if ~CtrlVar.AGlenisElementBased
            AGlenint=AGlen*fun;
            nint=n*fun;
        else
            AGlenint=AGlen;
            nint=n;
        end
        
        D=2*AGlenint.*(rhoint.*g).^nint./(nint+2);
        
        ds0dx=zeros(MUA.Nele,1); ds0dy=zeros(MUA.Nele,1);
        ds1dx=zeros(MUA.Nele,1); ds1dy=zeros(MUA.Nele,1);
        
        
        
        % derivatives for all elements at this integration point
        for I=1:MUA.nod
            
            ds0dx=ds0dx+Deriv(:,1,I).*s0nod(:,I);
            ds0dy=ds0dy+Deriv(:,2,I).*s0nod(:,I);
            ds1dx=ds1dx+Deriv(:,1,I).*s1nod(:,I);
            ds1dy=ds1dy+Deriv(:,2,I).*s1nod(:,I);
            
        end
        
        gradSurf1=sqrt(abs(ds1dx.*ds1dx+ds1dy.*ds1dy))+eps;
        gradSurf0=sqrt(abs(ds0dx.*ds0dx+ds0dy.*ds0dy))+eps;
        
        detJw=detJ*weights(Iint);
        
        for I=1:MUA.nod
            if ~OnlyR
                for J=1:MUA.nod
                    
                    deltahterm=fun(I).*fun(J);
                    
                    lf1=dt*theta*(n+2)*D.*(gradSurf1.^(n-1)).*(h1int.^(n+1))...
                        .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)).*fun(J);
                    
                    temp=dt*theta*D.*h1int.^(n+2).*(Deriv(:,1,I).*Deriv(:,1,J)+Deriv(:,2,I).*Deriv(:,2,J));
                    
                    lf2=(gradSurf1.^(n-1)).* temp;
                    lf3=(n-1)*(gradSurf1.^(n-3)).* temp.*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)) ;
                    
                    d1d1(:,I,J)=d1d1(:,I,J)+(deltahterm+lf1+lf2+lf3).*detJw;
                    
                end
            end
            dhterm=(h1int-h0int).*fun(I);
            
            q0term=dt*(1-theta)*D.*gradSurf0.^(n-1).*h0int.^(n+2).*(ds0dx.*Deriv(:,1,I)+ds0dy.*Deriv(:,2,I));
            q1term=dt*theta*    D.*gradSurf1.^(n-1).*h1int.^(n+2).*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I));
            
            a0term=-dt*(1-theta)*a0int.*fun(I);
            a1term=-dt*theta*a1int.*fun(I);
            
            % R=T-F
            % K du = -R
            %b1(:,I)=b1(:,I)+(dhterm+a0term+a1term+q0term+q1term).*detJw;
            t1(:,I)=t1(:,I)+ dhterm.*detJw;
            f1(:,I)=f1(:,I)-(a0term+a1term+q0term+q1term).*detJw;
            
        end
        %%-]
    end
    
end
% assemble right-hand side

T=sparseUA(neq,1);
for I=1:MUA.nod
    T=T+sparseUA(MUA.connectivity(:,I),ones(MUA.Nele,1),t1(:,I),neq,1);
end

F=sparseUA(neq,1);
for I=1:MUA.nod
    F=F+sparseUA(MUA.connectivity(:,I),ones(MUA.Nele,1),f1(:,I),neq,1);
end

R=T-F  ;

if ~OnlyR
    for I=1:MUA.nod
        for J=1:MUA.nod
            K=K+sparseUA(MUA.connectivity(:,I),MUA.connectivity(:,J),d1d1(:,I,J),neq,neq);
        end
    end
else
    K=[];
end

%K=(K+K')/2;

if CtrlVar.InfoLevelCPU>=10 ;  fprintf(CtrlVar.fidlog,' MatrixAssemblySSHEETtransient2HD in %-g sec. \n',toc(tCPU)) ; end



end