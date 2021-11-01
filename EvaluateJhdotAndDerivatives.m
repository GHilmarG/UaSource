function [Jhdot,duJhdot,dvJhdot,dhJhdot]=EvaluateJhdotAndDerivatives(UserVar,CtrlVar,MUA,F,BCs,Meas)


ndim=2; dof=1; neq=dof*MUA.Nnodes;

anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

dhdtMeasnod=reshape(Meas.dhdt(MUA.connectivity,1),MUA.Nele,MUA.nod);


% [points,weights]=sample('triangle',MUA.nip,ndim);

JhdotIntSum=zeros(MUA.Nele,1);
duJhdotIntSum=zeros(MUA.Nele,MUA.nod);
dvJhdotIntSum=zeros(MUA.Nele,MUA.nod);
dhJhdotIntSum=zeros(MUA.Nele,MUA.nod);

Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);
dhdtErr=sqrt(spdiags(Meas.dhdtCov));
dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),MUA.Nele,MUA.nod);


% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    aint=anod*fun;
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    dhdtMeasint=dhdtMeasnod*fun;
    dhdtErrInt=dhdtErrnod*fun;
    
    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);
    dudx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod
        
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
        
    end
    
    detJw=detJ*MUA.weights(Iint);
    
    JhdotIntSum=JhdotIntSum+((aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy)-dhdtMeasint)./dhdtErrInt).^2 .*detJw/2/Area; % element variable
    
    for Inod=1:MUA.nod
        
        
        
        duJhdotInt=-((aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy)-dhdtMeasint)./dhdtErrInt)...
            .*((dhdx.*fun(Inod)+hint.*Deriv(:,1,Inod))./dhdtErrInt)...
            .*detJw/Area;
        
        
        dvJhdotInt=-((aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy)-dhdtMeasint)./dhdtErrInt)...
            .*((dhdy.*fun(Inod)+hint.*Deriv(:,2,Inod))./dhdtErrInt)...
            .*detJw/Area;
        
        dhJhdotInt=-((aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy)-dhdtMeasint)./dhdtErrInt)...
            .*((dudx.*fun(Inod)+uint.*Deriv(:,1,Inod)+dvdy.*fun(Inod)+vint.*Deriv(:,2,Inod))./dhdtErrInt)...
            .*detJw/Area;
        
        
        duJhdotIntSum(:,Inod)=duJhdotIntSum(:,Inod)+duJhdotInt;
        dvJhdotIntSum(:,Inod)=dvJhdotIntSum(:,Inod)+dvJhdotInt;
        dhJhdotIntSum(:,Inod)=dhJhdotIntSum(:,Inod)+dhJhdotInt;
        
        
    end
end

% assemble right-hand side

Jhdot=sum(JhdotIntSum) ;

duJhdot=sparseUA(neq,1);
dvJhdot=sparseUA(neq,1);
dhJhdot=sparseUA(neq,1);
for Inod=1:MUA.nod
    
    duJhdot=duJhdot+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),duJhdotIntSum(:,Inod),neq,1);
    dvJhdot=dvJhdot+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),dvJhdotIntSum(:,Inod),neq,1);
    dhJhdot=dhJhdot+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),dhJhdotIntSum(:,Inod),neq,1);
    
end


%[hL,hRhs]=createLh(MUA.Nnodes,BCs.dhdtFixedNode,BCs.dhdtFixedValue,BCs.dhdtTiedNodeA,BCs.dhdtTiedNodeB);
%CtrlVar.SymmSolver='AugmentedLagrangian';
%x0=zeros(MUA.Nnodes,1) ; y0=hRhs*0;

%[dhdt,dhdtlambda]=solveKApeSymmetric(MUA.M,hL,Jhdot,hRhs,x0,y0,CtrlVar);

%[duJhdot,dhJhdot]=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,duJhdot,dhJhdot);


duJhdot=full(duJhdot);
dvJhdot=full(dvJhdot);
dhJhdot=full(dhJhdot);


end