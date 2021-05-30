function [t1,f1,d1d1]=SSHEETintAssembly(Iint,CtrlVar,MUA,AGlen,n,C,m,rhonod,g,s0nod,h0nod,s1nod,h1nod,a0nod,a1nod,dt,OnlyR)



% Notes:
%
%   -The natural BCs are no-flux
%   -Here the density has not (yet) been included in the mass conservation equation. (Doing so would be very simple, just multiply all terms by rho)
%

narginchk(17,17)
ndim=2;
theta=CtrlVar.theta;
%    [points,weights]=sample('triangle',MUA.nip,ndim);

if ~OnlyR
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
else
    d1d1=[];
end

t1=zeros(MUA.Nele,MUA.nod); f1=zeros(MUA.Nele,MUA.nod);

%% [--
fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % MUA.nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points

Deriv=MUA.Deriv(:,:,:,Iint);
detJ=MUA.DetJ(:,Iint);
%     % Deriv : MUA.Nele x dof x MUA.nod
%  detJ : MUA.Nele

% values at integration point


h0int=h0nod*fun;
h1int=h1nod*fun;
a0int=a0nod*fun;
a1int=a1nod*fun;
rhoint=rhonod*fun;


AGlenint=AGlen*fun;
nint=n*fun;

Cint=C*fun;
mint=m*fun;

Dd=2*AGlenint.*(rhoint.*g).^nint./(nint+2);
Db=Cint.*(rhoint.*g).^mint;


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

detJw=detJ*MUA.weights(Iint);

for I=1:MUA.nod
    if ~OnlyR
        for J=1:MUA.nod
            
            
            % the weight function is index I, the perturbation in h is index I
            % so replace \Delta h with fun(I) and \p_x \Delta h with Deriv(:,1,J)
            
            deltahterm=fun(I).*fun(J);
            
            lf1d=dt*theta*Dd.*(nint+2).*(gradSurf1.^(nint-1)).*(h1int.^(nint+1))...
                .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)).*fun(J);    % (2)
            lf2d=dt*theta*Dd.*h1int.^(nint+2).*(gradSurf1.^(nint-1)).*(Deriv(:,1,I).*Deriv(:,1,J)+Deriv(:,2,I).*Deriv(:,2,J));    % (1)
            lf3d=dt*theta.*Dd.*(nint-1).*h1int.^(nint+2).* gradSurf1.^(nint-3) .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)).*(ds1dx.*Deriv(:,1,J)+ds1dy.*Deriv(:,2,J)) ;  % (3)
            
            
            
            lf1b=dt*theta*Db.*(mint+1).*(gradSurf1.^(mint-1)).*(h1int.^mint)...
                .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)).*fun(J);   % (2)
            lf2b=dt*theta*Db.*h1int.^(mint+1).*(gradSurf1.^(mint-1)).*(Deriv(:,1,I).*Deriv(:,1,J)+Deriv(:,2,I).*Deriv(:,2,J));    % (1)
            lf3b=dt*theta*Db.*(mint-1).*h1int.^(mint+1).* gradSurf1.^(mint-3) .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)).*(ds1dx.*Deriv(:,1,J)+ds1dy.*Deriv(:,2,J)) ;  % (3)
            
            
            
            d1d1(:,I,J)=d1d1(:,I,J)+(deltahterm+lf1d+lf2d+lf3d+lf1b+lf2b+lf3b).*detJw;
            
        end
    end
    
    dhterm=(h1int-h0int).*fun(I);
    
    qd0term=dt*(1-theta)*Dd.* (gradSurf0.^(nint-1)) .* (h0int.^(nint+2)) .*(ds0dx.*Deriv(:,1,I)+ds0dy.*Deriv(:,2,I));
    qd1term=dt*theta*    Dd.* (gradSurf1.^(nint-1)) .* (h1int.^(nint+2)) .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I));
    
    qb0term=dt*(1-theta)*Db.* (gradSurf0.^(mint-1)) .* (h0int.^(mint+1)) .*(ds0dx.*Deriv(:,1,I)+ds0dy.*Deriv(:,2,I));
    qb1term=dt*theta*    Db.* (gradSurf1.^(mint-1)) .* (h1int.^(mint+1)) .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I));
    
    a0term=-dt*(1-theta)*a0int.*fun(I);
    a1term=-dt*theta*a1int.*fun(I);
    
    % R=T-F
    % K du = -R
    %b1(:,I)=b1(:,I)+(dhterm+a0term+a1term+q0term+q1term).*detJw;
    t1(:,I)=t1(:,I)+ dhterm.*detJw;
    f1(:,I)=f1(:,I)-(a0term+a1term+qd0term+qd1term+qb0term+qb1term).*detJw;
    
end

%%-]


end