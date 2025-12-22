


function [Jhdot,duJhdot,dvJhdot,dhJhdot]=EvaluateJhdotAndDerivatives(UserVar,CtrlVar,MUA,F,BCs,Meas)

%%  Provides cost function and derivatives with respect to $u$, $v$, and $h$, of the cost function term involving $\dot{h}$
%
% $$J_{\dot{h}} = \| \dot{h} - \hat{\dot{h}} \| $$
%
% where
%
% $$\dot{h} = a - ( \partial_x (u h ) + \partial (v h) ) $$ 
%
% and
%
% $$ J_{\dot{h}} = \frac{1}{2 \mathcal{A}} \int \! \int \left  (\dot{h} - \hat{\dot{h}} \right )^2 \; dx \, dy $$ 
%
%
% $$ d_u J_{\dot{h}} = \frac{1}{\mathcal{A}} \int \! \int (\dot{h} - \hat{\dot{h}} ) \, (\partial_u \dot{h} ) \; \delta u \; dx \, dy $$ 
%
%
% $$ \partial_u \dot{h} = \partial_ u  (  a - ( \partial_x (u h ) + \partial (v h) )) = - \partial_ u  \,  \partial_x (u h ) $$ 
%
% or
%
% $$ \partial_u \dot{h} = - \partial_x  \,  \partial_u (u h ) = -\partial_x (\delta u \, h) = - h \, \partial_x \delta u - \delta u \, \partial_x h $$ 
%
% $$ d_u J_{\dot{h}} = -\frac{1}{\mathcal{A}} \int \! \int (\dot{h} - \hat{\dot{h}} ) \, ( h \partial_x \delta u - \delta u \, \partial_x h ) \; dx \, dy $$ 
%
% Note: if F.dhdt is available, then one can simply calculate 
%
% $$ J_{\dot{h}} =  \frac{1}{2 \mathcal{A}}  (\dot{h}-\hat{\dot{h}})' \, M \, (\dot{h}-\hat{\dot{h}}) $$
%
% where $\dot{h}$=F.dhdt. Doing so ensures that $h$ boundary conditions have been accounted for.
%
% 
%
%
%%


ndim=2; dof=1; neq=dof*MUA.Nnodes;

anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

[~,F.dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs) ; 

if ~isempty(F.dhdt) || ~isnan(F.dhdt)
    dhdtnod=reshape(F.dhdt(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    dhdtnod=nan;
end

dhdtMeasnod=reshape(Meas.dhdt(MUA.connectivity,1),MUA.Nele,MUA.nod);


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


     if ~isnan(dhdtnod)
        hdot=dhdtnod*fun;
    else
        hdot=aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy) ;
     end


    
    detJw=detJ*MUA.weights(Iint);
    
  
    JhdotIntSum=JhdotIntSum+((hdot-dhdtMeasint)./dhdtErrInt).^2 .*detJw/2/Area; 
    
    for Inod=1:MUA.nod
        
        % hdot=aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy) ; 
        
        duJhdotInt=-((hdot-dhdtMeasint)./dhdtErrInt)...
            .*((dhdx.*fun(Inod)+hint.*Deriv(:,1,Inod))./dhdtErrInt)...
            .*detJw/Area;
        
        
        dvJhdotInt=-((hdot-dhdtMeasint)./dhdtErrInt)...
            .*((dhdy.*fun(Inod)+hint.*Deriv(:,2,Inod))./dhdtErrInt)...
            .*detJw/Area;
        
        dhJhdotInt=-((hdot-dhdtMeasint)./dhdtErrInt)...
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


duJhdot=full(duJhdot);
dvJhdot=full(dvJhdot);
dhJhdot=full(dhJhdot);


%% If F.dhdt is available this should give the same answer
% dhdtErr=sqrt(spdiags(Meas.dhdtCov)) ;  dhdtres=(F.dhdt-Meas.dhdt)./dhdtErr ;  JhdotTest=full(dhdtres'*MUA.M*dhdtres)/2/Area;
%%

end