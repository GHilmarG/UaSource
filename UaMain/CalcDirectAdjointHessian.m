


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

dudA=[]; dvdA=[];
dudB=[]; dvdB=[];
dudC=[]; dvdC=[];

%CtrlVar.Calculate.Geometry="bh-FROM-sBS" ; 

if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    tA=tic;
    [dudA,dvdA]=duvdAFunc(CtrlVar,MUA,F,BCs) ;  % this has been tested against finite-differences and is good
    tA=toc(tA);
    fprintf("dudA and dvdA sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    tB=tic;
    [dudB,dvdB]=duvdBFunc(CtrlVar,MUA,F,BCs) ;  % this has been tested against finite-differences and is good
    tB=toc(tB);
    fprintf("dudB and dvdB sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tB)
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    tC=tic;
    [dudC,dvdC]=duvdCFunc(CtrlVar,MUA,F,BCs) ; % this has been tested against finite-differences and is good
    tC=toc(tC);
    fprintf("dudC and dvdC sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tC)
end




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
if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)

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
xi=[dudA dudB dudC ; dvdA dvdB dvdC] ; % 2 Nnodes \times nP Nnodes where nP is the number of fields inverted for, e.g. 2 if inverting for A and C, 1 in only inverting for B

d2Jdqq=blkdiag(d2Iduu,d2Idvv);  % d2Iduv and d2Idvu are zeros, 2 Nnodes \times 2 Nnodes


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

    % Fxi=FindOrCreateFigure("xi=dq/dp") ; contourf(rot90(xi),LineStyle="none") ; axis equal ; colorbar
    % title("abs(xi)")
    % CM=cmocean('balanced',25,'pivot',0) ; colormap(Fxi,CM);
    % 
    % %set(gca,'ColorScale','log')
    % 
    % FindOrCreateFigure("d2Idpp") ; contourf(abs(rot90(d2Idpp)),LineStyle="none") ; axis equal ; colorbar
    % title("abs(d2Idpp)")
    % set(gca,'ColorScale','log')

    NodeTest=804;

    %% A 
    if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
        % du/dA  But I need log10 sensitivities
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

        % dv/dA
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

    end


      %% B
    if contains(CtrlVar.Inverse.InvertFor,"-B-")
        % du/dB  
          
        CtrlVar.Calculate.Geometry="bh-FROM-sBS" ; 

        dudB=dudB(:,NodeTest);
        dvdB=dvdB(:,NodeTest);
     
        B0=F.B;
        dB=1;

        Bp=B0;
        Bp(NodeTest)=Bp(NodeTest)+dB;
        F.B=Bp;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        up=F.ub; vp=F.vb;

        
        Bm=B0;
        Bm(NodeTest)=Bm(NodeTest)-dB;
        F.B=Bm;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        um=F.ub; vm=F.vb;

        F.B=B0;

        dudBpert=(up-um)/(2*dB) ;  dvdBpert=(vp-vm)/(2*dB) ;

        % dv/dB
        figBu=FindOrCreateFigure("du/dB comparision");
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudB,CreateNewFigure=false)  ; title("$du/dB$ sensitvity",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudBpert,CreateNewFigure=false) ; title("$du/dB$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,dudB-dudBpert,CreateNewFigure=false) ; title("$du/dB$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";

        figBv=FindOrCreateFigure("dv/dB comparision");
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdB,CreateNewFigure=false)  ; title("$dv/dB$ sensitvity",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdBpert,CreateNewFigure=false) ; title("$dv/dB$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,dvdB-dvdBpert,CreateNewFigure=false) ; title("$dv/dB$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";


        figBgrad=FindOrCreateFigure("du/dB gradient test") ;  clf(figBgrad)
        plot(dudB,dudBpert,"or") ;
        hold on
        axis equal
        AX=axis;
        plot([min(dudB) max(dudB)],[min(dudB) max(dudB)],"--k") ;
        axis equal tight ;
        xlabel(" $du/dB$",Interpreter="latex")  ;
        ylabel("Finite difference $du/dB$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')

        figBgrad=FindOrCreateFigure("dv/dB gradient test") ;  clf(figBgrad)
        plot(dvdB,dvdBpert,"or") ;
        hold on
        axis equal
        plot([min(dvdB) max(dvdB)],[min(dvdB) max(dvdB)],"--k") ;
        axis equal tight ;
        xlabel(" $dv/dB$",Interpreter="latex")  ;
        ylabel("Finite difference $dv/dB$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')





    end

    %% C
    if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
         % du/dC
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

        % dv/dC
        figBu=FindOrCreateFigure("du/dC comparision");
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


end