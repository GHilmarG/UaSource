

function [s,b,u,x]=AnalyticalOneDimentionalIceShelf(CtrlVar,MUA,F,hgl,ugl,xgl)
    %
    %
    % a=0.3;
    % n=3;
    % rho=910;
    % rhow=1030;
    % g=9.81/1000;
    % A=AGlenVersusTemp(-10);
    % hgl=2000;
    % ugl=300;
    
    
    a=F.as+F.ab ;
    n=F.n;
    rho=F.rho;
    rhow=F.rhow;
    g=F.g;
    A=F.AGlen;
    
    qgl=hgl*ugl;
    gamm=A.*(rho.*(1-rho./rhow).*g).^n./(4.^n);
    x=sort(MUA.coordinates(:,1)-xgl) ;
    
    
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