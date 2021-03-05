

function [s,b,u,x]=AnalyticalOneDimentionalIceShelf(CtrlVar,MUA,hgl,ugl,xgl,x,A,n,rho,rhow,a,g)
    %
    
    % sorts the outputs with x if x is requied as an output.
    % otherwise, s, b, u are ordered as the nodes in the mesh
    %
    %
    %
    % Examples: 
    %
    %
    % Use default values, and x values based on MUA.coordinates(:,1)
    %    [s,b,u,x]=AnalyticalOneDimentionalIceShelf([],MUA)
    %
    % Solution at x using defaults values
    %    [s,b,u]=AnalyticalOneDimentionalIceShelf([],x)
    %
    %%
    
    if ~ ( (nargin == 2) || (nargin==13) )
        error('AnalyticalOneDimentionalIceShelf:IncorrectNumberOfInputArguments','Incorrect number of input arguments')
    end
    
    if nargin==2
        % fprintf('AnalyticalOneDimentionalIceShelf using default values.\n')
        a=0.3;
        n=3;
        rho=910;
        rhow=1030;
        g=9.81/1000;
        A=AGlenVersusTemp(-10);
        hgl=1000;
        ugl=300;
        xgl=0 ;
        if isstruct(MUA)
            x=MUA.coordinates(:,1) ;
        else
            x=MUA;
        end
        
    end

    
    
    x=x-xgl ; 
    
    qgl=hgl*ugl;
    gamm=A.*(rho.*(1-rho./rhow).*g).^n./(4.^n);

    
    K=qgl.^(n + 1).*(a./hgl.^(n + 1)-gamm);
    
    if a==0
        h=(hgl.^(-(n+1))+gamm.*(n+1).*x./qgl).^(-1/(n+1));
    else
        h=( (gamm+K./((qgl+a.*x).^(n+1)))./a).^(-1./(n+1));
    end
    
    u=(a.*x+qgl)./h;
    s=(1-rho./rhow).*h;
    b=s-h;
    
    x=x+xgl ; 
    
    if nargout> 3
       
        [x,I]=sort(x);
        b=b(I) ;
        s=s(I);
        u=u(I); 
        
    end
    
    
end