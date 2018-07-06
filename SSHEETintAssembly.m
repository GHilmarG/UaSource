function [t1,f1,d1d1]=SSHEETintAssembly(Iint,CtrlVar,MUA,AGlen,n,rhonod,g,s0nod,h0nod,s1nod,h1nod,a0nod,a1nod,dt,OnlyR)
        
    ndim=2; 
    theta=CtrlVar.theta;
    [points,weights]=sample('triangle',MUA.nip,ndim);
    
    if ~OnlyR
        d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
    else
        d1d1=[];  
    end
    
    t1=zeros(MUA.Nele,MUA.nod); f1=zeros(MUA.Nele,MUA.nod);
    
    %% [--
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % MUA.nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        
        if isfield(MUA,'Deriv') && isfield(MUA,'DetJ')
            Deriv=MUA.Deriv(:,:,:,Iint);
            detJ=MUA.DetJ(:,Iint);
        else
            [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
        end
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
                    
                    lf1=dt*theta*(nint+2).*D.*(gradSurf1.^(nint-1)).*(h1int.^(nint+1))...
                        .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)).*fun(J);
                    
                    temp=dt*theta*D.*h1int.^(nint+2).*(Deriv(:,1,I).*Deriv(:,1,J)+Deriv(:,2,I).*Deriv(:,2,J));
                    
                    lf2=(gradSurf1.^(nint-1)).* temp;
                    lf3=(nint-1).*(gradSurf1.^(nint-3)).* temp.*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I)) ; 

                    %             lf2=dt*theta*D.*(gradSurf1.^(nint-1)).*(h1int.^(nint+2))...
                    %                 .*(Deriv(:,1,I).*Deriv(:,1,J)+Deriv(:,2,I).*Deriv(:,2,J));
                    %
                    %             lf3=dt*theta*(nint-1)*D.*(gradSurf1.^(nint-3)).*(h1int.^(nint+2))...
                    %                 .*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I))...
                    %                 .*(Deriv(:,1,I).*Deriv(:,1,J)+Deriv(:,2,I).*Deriv(:,2,J));
                    
                    d1d1(:,I,J)=d1d1(:,I,J)+(deltahterm+lf1+lf2+lf3).*detJw;
                    
                end
            end
            dhterm=(h1int-h0int).*fun(I);
            
            q0term=dt*(1-theta)*D.*gradSurf0.^(nint-1).*h0int.^(nint+2).*(ds0dx.*Deriv(:,1,I)+ds0dy.*Deriv(:,2,I));
            q1term=dt*theta*    D.*gradSurf1.^(nint-1).*h1int.^(nint+2).*(ds1dx.*Deriv(:,1,I)+ds1dy.*Deriv(:,2,I));
            
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