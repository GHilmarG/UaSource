function [wSurfInt,wBedInt]=calcVerticalSurfaceVelocity2(h,b,u,v,exx,eyy,xint,yint,coordinates,connectivity,nip)
    
    % calculates vertical velocity at both integration points and nodal points
    
    % NOTE: currently assumes vertical velocity zero at bed (does not work for ice shelfs)
    
    % get h at integration points
    
    [Nele,nod]=size(connectivity); ndim=2;
    [points,weights]=sample('triangle',nip,ndim);
    
    hnod=reshape(h(connectivity,1),Nele,nod);
    bnod=reshape(b(connectivity,1),Nele,nod);
    unod=reshape(u(connectivity,1),Nele,nod);
    vnod=reshape(v(connectivity,1),Nele,nod);
    
   % coox=reshape(coordinates(connectivity,1),Nele,nod);
   % cooy=reshape(coordinates(connectivity,2),Nele,nod);
    
    
   wBedInt=zeros(Nele,nip); wSurfInt=zeros(Nele,nip);
   
   
       
    for Iint=1:nip                           % loop over integration points
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        hint=hnod*fun;
        %bint(:,Iint)=bnod*fun;
        uint=unod*fun;
        vint=vnod*fun;
        
        [Deriv]=derivVector(coordinates,connectivity,nip,Iint); % Deriv=zeros(Nele,dof,nod);
        % deriv(1,:)=Deriv(Iele,nip,1,:)
        
        dbdxint=zeros(Nele,1) ; dbdyint=zeros(Nele,1) ;
        for I=1:nod
            dbdxint=dbdxint+Deriv(:,1,I).*bnod(:,I);
            dbdyint=dbdyint+Deriv(:,2,I).*bnod(:,I);
        end
         wBedInt(:,Iint)=uint.*dbdxint+vint.*dbdyint;
         wSurfInt(:,Iint)=-hint.*(exx(:,Iint)+eyy(:,Iint))+wBedInt(:,Iint);
    
    end
   
    
    % calculate w(b) for grounded ice
    % w(b)=u db/dx + v db/dy
    % must calculate u and v on integra
    % tion points and bed slopes
    % bed slopes
    
    
    
end
