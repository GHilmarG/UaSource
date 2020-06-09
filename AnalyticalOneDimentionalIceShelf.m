

function [s,b,u,x]=AnalyticalOneDimentionalIceShelf(CtrlVar,MUA,F,hgl,ugl,xgl,x,A,n,rho,rhow,a,g)
    %
    %
    % a=0.3;
    % n=3;
    % rho=910;
    % rhow=1030;
    % g=9.81/1000;
    % A=AGlenVersusTemp(-10);
    % hgl=1000;
    % ugl=300;
    
    if nargin<9
        a=F.as(1)+F.ab(1) ;
        n=F.n(1);
        rho=F.rho(1);
        rhow=F.rhow(1);
        g=F.g(1);
        A=F.AGlen(1);
    end
    
    qgl=hgl*ugl;
    gamm=A.*(rho.*(1-rho./rhow).*g).^n./(4.^n);
    
    if nargin<7 || isempty(x)
        x=sort(MUA.coordinates(:,1)-xgl) ;
    end
    
    K=qgl.^(n + 1).*(a./hgl.^(n + 1)-gamm);
    
    if a==0
        h=(hgl.^(-(n+1))+gamm.*(n+1).*x./qgl).^(-1/(n+1));
    else
        h=( (gamm+K./((qgl+a.*x).^(n+1)))./a).^(-1./(n+1));
    end
    
    u=(a.*x+qgl)./h;
    s=(1-rho./rhow).*h;
    b=s-h;
    
    
end