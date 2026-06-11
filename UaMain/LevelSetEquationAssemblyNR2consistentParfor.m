   

function [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistentParfor(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)

 
    narginchk(15,53)
    
    
    ndim=2; dof=1; neq=dof*MUA.Nnodes;
    
    theta=CtrlVar.LevelSetTheta;
    dt=CtrlVar.dt;
    CtrlVar.Tracer.SUPG.tau=CtrlVar.LevelSetSUPGtau;

    isL=CtrlVar.LSF.L ; isP=CtrlVar.LSF.P ; isT=CtrlVar.LSF.T ;

    f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    f1nod=reshape(f1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    
    c0nod=reshape(c0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    c1nod=reshape(c1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    qx0nod=reshape(qx0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    qy0nod=reshape(qy0(MUA.connectivity,1),MUA.Nele,MUA.nod);

    qx1nod=reshape(qx1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    qy1nod=reshape(qy1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
    %b1=zeros(MUA.Nele,MUA.nod);
    TG=zeros(MUA.Nele,MUA.nod);
    PG=zeros(MUA.Nele,MUA.nod);
    LG=zeros(MUA.Nele,MUA.nod);
    R=zeros(MUA.Nele,MUA.nod);
    RSUPG=zeros(MUA.Nele,MUA.nod);
    
    qx=zeros(MUA.Nele,MUA.nod);
    qy=zeros(MUA.Nele,MUA.nod);
    
    
    if CtrlVar.LevelSetSolutionMethod=="Newton Raphson"
        NR=1;
    else
        NR=0;
    end
    
    PGsum=PG; LGsum=LG ; TGsum=TG ; Rsum=R ; RSUPGsum=RSUPG ; qxSum=qx ; qySum=qy ; d1d1Sum=d1d1;
    Deriv=MUA.Deriv;
    detJ=MUA.DetJ;
    weights=MUA.weights;
    Nele=MUA.Nele;
    nod=MUA.nod;
    points=MUA.points;
    EleAreas=MUA.EleAreas;    
    
    
    parfor Iint=1:MUA.nip  %Integration points
        
        %Deriv=MUA.Deriv(:,:,:,Iint);
        %detJ=MUA.DetJ(:,Iint);
        % important to reduce memory transfer, get rid of MUA
        [PG,LG,TG,R,RSUPG,qx,qy,d1d1]=AssemblyLSFintLoop(Iint,f0nod,f1nod,u0nod,u1nod,v0nod,v1nod,qx0nod,qy0nod,qx1nod,qy1nod,c0nod,c1nod,CtrlVar,NR,theta,dt,ndim,Deriv(:,:,:,Iint),detJ(:,Iint),weights(Iint),Nele,nod,points,EleAreas);
        PGsum=PGsum+PG;
        LGsum=LGsum+LG;
        TGsum=TGsum+TG;
        Rsum=Rsum+R;
        RSUPGsum=RSUPGsum+RSUPG;
        qxSum=qxSum+qx;
        qySum=qySum+qy;
        d1d1Sum=d1d1Sum+d1d1;
        
        
    end
    
    PG=PGsum; LG=LGsum ; TG=TGsum;  R=Rsum ; RSUPG=RSUPGsum ; qx=qxSum ; qy=qySum ; d1d1=d1d1Sum ;
    
    
    
    Pv=sparseUA(neq,1);
    Lv=sparseUA(neq,1);
    Tv=sparseUA(neq,1);
    Qx=sparseUA(neq,1);
    Qy=sparseUA(neq,1);
    Rv=sparseUA(neq,1);
    RSUPGv=sparseUA(neq,1);
    
    for Inod=1:MUA.nod
        Pv=Pv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),PG(:,Inod),neq,1);
        Lv=Lv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),LG(:,Inod),neq,1);
        Tv=Tv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),TG(:,Inod),neq,1);
        Rv=Rv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),R(:,Inod),neq,1);
        RSUPGv=RSUPGv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),RSUPG(:,Inod),neq,1);
        Qx=Qx+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qx(:,Inod),neq,1);
        Qy=Qy+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qy(:,Inod),neq,1);
    end
    
    rh=isL*Lv+isP*Pv+isT*Tv+RSUPGv; 
    
    if nargout>2
        Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
        istak=0;
        
        for Inod=1:MUA.nod
            for Jnod=1:MUA.nod
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
                Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
                Xval(istak+1:istak+MUA.Nele)=d1d1(:,Inod,Jnod);
                istak=istak+MUA.Nele;
            end
        end
        
        kv=sparseUA(Iind,Jind,Xval,neq,neq);
    end
    
end