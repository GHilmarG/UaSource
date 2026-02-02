


function d2Idpp = CalcDirectAdjointHessian(UserVar,CtrlVar,RunInfo,MUA,F,BCs,l,BCsAdjoint,d2Iduu,d2Idvv,d2Idhdothdot,uAdjoint,vAdjoint)

narginchk(13,13)


%% Calculates some of the terms of the Hessian using the direct-adjoint approach.
%
% Here this is done for the misfit term and later the Hessian of the regularization term is added (this is easy).
%
% The tricky part is to do this for the misfit term.
%
%
% $$
%   H_{ij} = \frac{\partial^2 J}{\partial p_i \, \partial p_j}
%   + \Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial p_j}
%   +\frac{\partial^2 J}{\partial q_k\, \partial q_m} \xi_{ki} \, \xi_{mj}
%   +\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial q_m} \xi_{ki} \, \xi_{mj}
%   +\frac{\partial^2 J}{\partial p_i \, \partial q_k} \xi_{kj}
%   +\Psi_n \frac{\partial^2 F_n}{\partial p_i \, \partial q_k} \, \xi_{kj}
%   +\frac{\partial^2 J}{\partial q_k \, \partial p_j} \xi_{ki}
%   +\Psi_n \frac{\partial^2 F_n}{\partial q_k \, \partial p_j} \xi_{ki}
% $$
%
%
% 
% Currently only for u, v as q variables
%
%
% $$2 I= u M_u u + v M_v v + \dot{h} M_{\dot{h}} \dot{h} $$
%
% $$ p =\left ( \begin{array}{c}  A \\ B \\ C  \end{array} \right ) $$
%
% $$ q =\left ( \begin{array}{c}  u \\ v \\ \dot{h}  \end{array} \right ) $$
%
% $$ \xi = \frac{\partial q}{\partial p} = \left ( \begin{array}{cc}
%  \frac{\partial u}{ \partial A} & \frac{\partial u}{\partial C}  \\
%  \frac{\partial v}{\partial A } & \frac{\partial v}{\partial C}  \\
%  \frac{\partial \dot{h}}{\partial A} & \frac{\partial \dot{h}}{\partial C}  \\
% \end{array} \right ) $$
%
% $$ \frac{\partial^2 J}{\partial q\, \partial q} = \left ( \begin{array}{ccc}
%                                                    M_u & 0 & 0 \\
%                                                    0 & M_v & 0 \\
%                                                    0 & 0 & M_{\dot{h}}
% \end{array} \right ) $$
%
%
% If the forward model is
%
% $$ F(q(p),p) = 0 $$
%
% where $q$ are output variables and $p$ model parameters, then
%
% $$ \partial F/\partial q \; \partial q / \partial p + \partial F / \partial p = 0 $$
%
% or
%
% $$ \frac{\partial F}{\partial q} \; \frac{\partial q }{ \partial p} = - \frac{\partial F }{ \partial p}  $$
%
% which can be solved for the sensitives
%
% $$ \xi_{ij} : = \frac{\partial q_i}{\partial p_j} $$ 
% 
%%

%% Get the sensitivity matrices


tA=tic;
[dudA,dvdA]=duvdAFunc(CtrlVar,MUA,F,BCs) ;  % this has been tested against finite-differences and is good
tA=toc(tA);
fprintf("dudA and dvdA sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)

tC=tic;
[dudC,dvdC]=duvdCFunc(CtrlVar,MUA,F,BCs) ; % this has been tested against finite-differences and is good
tC=toc(tC);
fprintf("dudC and dvdC sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tC)


%%
% log10 sensitivities
%
% du/dA=du/dx  dx/dA
%
% x=ln(A)
%
% du/dA=du/dx  d(ln(A))/dA  = du/dx   1/A
%
% Therefore
%
% du/d(ln(A)) = A du/dA
%
% or
%
% du/d(ln(A)) = log(10) A du/dA
%
%
if contains(CtrlVar.Inverse.Regularize.Field,"logaglen",IgnoreCase=true)

    scale=log(10)*F.AGlen';  % this has to be a row vector
    dudA=dudA.*scale ; % using implicit expansion
    dvdA=dvdA.*scale ; % using implicit expansion
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
  
    scale=log(10)*F.C';  % this has to be a row vector
    dudC=dudC.*scale ; % using implicit expansion
    dvdC=dvdC.*scale ; % using implicit expansion

end



%%

xi=[dudA dudC ; dvdA dvdC];

d2Jdqq=blkdiag(d2Iduu,d2Idvv);  % d2Iduv and d2Idvu are zeros


d2Idpp=xi'*(d2Jdqq*xi);

d2Idpp=0.5*(d2Idpp+d2Idpp');




TestSensitivites=false;

if TestSensitivites

    %% Test

    S=d2Idpp;
    upper=max(abs(S(:)));
    lower=min(abs(S(:)));
    Threshold=upper/1e10;
    S(abs(S)<Threshold)=0;
    S=sparse(S);
    figure(10); spy(S) ; title("d2Idpp")

    Fxi=FindOrCreateFigure("xi=dq/dp") ; contourf(rot90(xi),LineStyle="none") ; axis equal ; colorbar
    title("abs(xi)")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(Fxi,CM);

    %set(gca,'ColorScale','log')

    FindOrCreateFigure("d2Idpp") ; contourf(abs(rot90(d2Idpp)),LineStyle="none") ; axis equal ; colorbar
    title("abs(d2Idpp)")
    set(gca,'ColorScale','log')

    NodeTest=100;

    %% du/dA  But I need log10 sensitivities
    dudA=dudA(:,NodeTest);
    dvdA=dvdA(:,NodeTest);


    A0=F.AGlen;
    dA=F.AGlen(NodeTest)*0.00001;

    Ap=A0;
    Ap(NodeTest)=Ap(NodeTest)+dA;
    F.AGlen=Ap;
    [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
    up=F.ub; vp=F.vb;

    F.AGlen=A0;
    Am=A0;
    Am(NodeTest)=Am(NodeTest)-dA;
    F.AGlen=Am;
    [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
    um=F.ub; vm=F.vb;

    F.AGlen=A0;

    dudApert=(up-um)/(2*dA) ;  dvdApert=(vp-vm)/(2*dA) ;

    %% dv/dA
    figAu=FindOrCreateFigure("du/dA comparision");
    T=tiledlayout("flow");

    T1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dudA,CreateNewFigure=false)  ; title("$du/dA$ sensitvity",Interpreter="latex") ; subtitle("")
    title(cbar,"")

    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dudApert,CreateNewFigure=false) ; title("$du/dA$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dudA-dudApert,CreateNewFigure=false) ; title("$du/dA$ differences",Interpreter="latex") ; subtitle("")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

    T.Padding="loose";   T.TileSpacing="tight";

    figAv=FindOrCreateFigure("dv/dA comparision");
    T=tiledlayout("flow");

    T1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dvdA,CreateNewFigure=false)  ; title("$dv/dA$ sensitvity",Interpreter="latex") ; subtitle("")
    title(cbar,"")

    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dvdApert,CreateNewFigure=false) ; title("$dv/dA$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dvdA-dvdApert,CreateNewFigure=false) ; title("$dv/dA$ differences",Interpreter="latex") ; subtitle("")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

    T.Padding="loose";   T.TileSpacing="tight";


    %% du/dC
    dudC=dudC(:,NodeTest);
    dvdC=dvdC(:,NodeTest);


    C0=F.C;
    dC=F.C(NodeTest)*0.00001;

    Cp=C0;
    Cp(NodeTest)=Cp(NodeTest)+dC;

    F.C=Cp;
    [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
    up=F.ub; vp=F.vb;

    F.C=C0;
    Cm=C0;
    Cm(NodeTest)=Cm(NodeTest)-dC;
    F.C=Cm;
    [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
    um=F.ub; vm=F.vb;

    F.C=C0;

    dudCpert=(up-um)/(2*dC) ;  dvdCpert=(vp-vm)/(2*dC) ;

    %% dv/dC
    figAu=FindOrCreateFigure("du/dC comparision");
    T=tiledlayout("flow");

    T1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dudC,CreateNewFigure=false)  ; title("$du/dC$ sensitvity",Interpreter="latex") ; subtitle("")
    title(cbar,"")

    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dudCpert,CreateNewFigure=false) ; title("$du/dC$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dudC-dudCpert,CreateNewFigure=false) ; title("$du/dC$ differences",Interpreter="latex") ; subtitle("")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

    T.Padding="loose";   T.TileSpacing="tight";

    figAv=FindOrCreateFigure("dv/dC comparision");
    T=tiledlayout("flow");

    T1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dvdC,CreateNewFigure=false)  ; title("$dv/dC$ sensitvity",Interpreter="latex") ; subtitle("")
    title(cbar,"")

    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dvdCpert,CreateNewFigure=false) ; title("$dv/dA$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

    T3=nexttile;
    UaPlots(CtrlVar,MUA,F,dvdC-dvdCpert,CreateNewFigure=false) ; title("$dv/dA$ differences",Interpreter="latex") ; subtitle("")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

    T.Padding="loose";   T.TileSpacing="tight";


end


end