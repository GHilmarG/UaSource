
function [dudA,dvdA,dhdA,dudB,dvdB,dhdB,dudC,dvdC,dhdC]=duv_hdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs)


narginchk(7,7)
nargoutchk(9,9)

dudA=[]; dvdA=[]; dhdA=[];
dudB=[]; dvdB=[]; dhdB=[];
dudC=[]; dvdC=[]; dhdC=[];

CtrlVar.uvAssembly.ZeroFields=false;
CtrlVar.uvMatrixAssembly.Ronly=false;
[~,KdFuvduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);


if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    tA=tic;
    [dudA,dvdA,dhdA]=duvhdotdAFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good, also for dhdotdA
    tA=toc(tA);
    fprintf("A sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    tB=tic;
    [dudB,dvdB,dhdB]=duvdBFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good
    tB=toc(tB);
    fprintf("B sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tB)
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    tC=tic;
    [dudC,dvdC,dhdC]=duvdCFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ; % this has been tested against finite-differences and is good
    tC=toc(tC);
    fprintf("C sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tC)
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
    dhdA=dhdA.*scale ;

end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)

    scale=log(10)*F.C';  % this has to be a row vector
    dudC=dudC.*scale ;   % using implicit expansion
    dvdC=dvdC.*scale ;   % using implicit expansion
    dhdC=dhdC.*scale ;

end

TestSensitivites=false;

if TestSensitivites

    %% Test


    Funperturbed=F;
    NodeTest=804;
    NodeTest=1200;
    %NodeTest=1500;

    %% A
    if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)


        % du/dA  : these are now log10 sensitivities
        dudA=dudA(:,NodeTest);
        dvdA=dvdA(:,NodeTest);

        if ~isempty(dhdA)
            dhdA=dhdA(:,NodeTest);
        end

        %%
        A0=F.AGlen;
        dA=F.AGlen(NodeTest)*0.01;

        Ap=A0;
        Ap(NodeTest)=Ap(NodeTest)+dA;
        F.AGlen=Ap;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
        up=F.ub; vp=F.vb; dhdtp=dhdt;

        F.AGlen=A0;
        Am=A0;
        Am(NodeTest)=Am(NodeTest)-dA;
        F.AGlen=Am;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
        um=F.ub; vm=F.vb; dhdtm=dhdt;

        F.AGlen=A0;

        dudApert=(up-um)/(2*dA) ;
        dvdApert=(vp-vm)/(2*dA) ;
        dhdotdApert=(dhdtp-dhdtm)/(2*dA) ;

        scale=log(10)*A0;
        dudApert=dudApert.*scale ;
        dvdApert=dvdApert.*scale ;
        dhdotdApert=dhdotdApert.*scale ;


        % dv/dA
        figAu=FindOrCreateFigure("du/dA comparision");
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudA,CreateNewFigure=false)  ; title("$du/dA$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudApert,CreateNewFigure=false) ; title("$du/dA$ finite differences ($\log$ scale) ",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,dudA-dudApert,CreateNewFigure=false) ; title("$du/dA$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";

        figAv=FindOrCreateFigure("dv/dA comparision");
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdA,CreateNewFigure=false)  ; title("$dv/dA$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdApert,CreateNewFigure=false) ; title("$dv/dA$ finite differences ($\log$ scale)",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,dvdA-dvdApert,CreateNewFigure=false) ; title("$dv/dA$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";


        % dhdot sensitivity to A

        if ~isempty(dhdA)


            % dhdot/dA
            figAhdot=FindOrCreateFigure("dhdot/dA comparision"); clf(figAhdot)
            T=tiledlayout("flow");

            T1=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdA,CreateNewFigure=false)  ; title("$d\dot{h}/dA$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
            title(cbar,"")

            T2=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdotdApert,CreateNewFigure=false) ; title("$d\dot{h}/dA$ finite differences ($\log$ scale) ",Interpreter="latex") ; subtitle("") ; title(cbar,"")

            T3=nexttile;
            UaPlots(CtrlVar,MUA,F,dhdA-dhdotdApert,CreateNewFigure=false) ; title("$d\dot{h}/dA$ differences",Interpreter="latex") ; subtitle("")
            CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

            T.Padding="loose";   T.TileSpacing="tight";


        end


        figAgradu=FindOrCreateFigure("du/dA gradient test") ;  clf(figAgradu)
        plot(dudA,dudApert,"or") ;
        hold on
        axis equal
        AX=axis;
        plot([min(dudA) max(dudA)],[min(dudA) max(dudA)],"--k") ;
        axis equal tight ;
        xlabel(" $du/dA$",Interpreter="latex")  ;
        ylabel("Finite difference $du/dA$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')

        figAgradv=FindOrCreateFigure("dv/dA gradient test") ;  clf(figAgradv)
        plot(dvdA,dvdApert,"or") ;
        hold on
        axis equal
        plot([min(dvdA) max(dvdA)],[min(dvdA) max(dvdA)],"--k") ;
        axis equal tight ;
        xlabel(" $dv/dA$",Interpreter="latex")  ;
        ylabel("Finite difference $dv/dA$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')


        if ~isempty(dhdA)

            figAgradh=FindOrCreateFigure("dhdot/dA gradient test") ;  clf(figAgradh)
            plot(dhdA,dhdotdApert,"or") ;
            hold on
            axis equal
            plot([min(dhdA) max(dhdA)],[min(dhdA) max(dhdA)],"--k") ;
            axis equal tight ;
            xlabel(" $d\dot{h}/dA$",Interpreter="latex")  ;
            ylabel("Finite difference $d\dot{h}/dA$",Interpreter="latex")
            ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
            axis on ; axis equal tight ; box off
            title("Comparision between adjoint and finite-differences gradient calculations")
            set(gcf,'Color','white')

        end




        %%
    end


    %% B
    if contains(CtrlVar.Inverse.InvertFor,"-B-")
        % du/dB
        F=Funperturbed;

        CtrlVar.Calculate.Geometry="bh-FROM-sBS" ;

        dudB=dudB(:,NodeTest);
        dvdB=dvdB(:,NodeTest);
        if ~isempty(dhdB)
            dhdB=dhdB(:,NodeTest);
        end

        B0=F.B;
        dB=1;

        Bp=B0;
        Bp(NodeTest)=Bp(NodeTest)+dB;
        F.B=Bp;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
        up=F.ub; vp=F.vb; dhdtp=dhdt;


        Bm=B0;
        Bm(NodeTest)=Bm(NodeTest)-dB;
        F.B=Bm;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
        um=F.ub; vm=F.vb; dhdtm=dhdt;

        F=Funperturbed;

        dudBpert=(up-um)/(2*dB) ;
        dvdBpert=(vp-vm)/(2*dB) ;
        dhdotdBpert=(dhdtp-dhdtm)/(2*dB) ;

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


        if ~isempty(dhdB)


            % dhdot/dB
            figBhdot=FindOrCreateFigure("dhdot/dB comparision"); clf(figBhdot)
            T=tiledlayout("flow");

            T1=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdB,CreateNewFigure=false)  ; title("$d\dot{h}/dB$ sensitvity ",Interpreter="latex") ; subtitle("")
            title(cbar,"")

            T2=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdotdBpert,CreateNewFigure=false) ; title("$d\dot{h}/dB$ finite differences ",Interpreter="latex") ; subtitle("") ; title(cbar,"")

            T3=nexttile;
            UaPlots(CtrlVar,MUA,F,dhdB-dhdotdBpert,CreateNewFigure=false) ; title("$d\dot{h}/dB$ differences",Interpreter="latex") ; subtitle("")
            CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

            T.Padding="loose";   T.TileSpacing="tight";


        end


        figBgradu=FindOrCreateFigure("du/dB gradient test") ;  clf(figBgradu)
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

        figBgradv=FindOrCreateFigure("dv/dB gradient test") ;  clf(figBgradv)
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


        if ~isempty(dhdB)

            figBgradh=FindOrCreateFigure("dhdot/dB gradient test") ;  clf(figBgradh)
            plot(dhdB,dhdotdBpert,"or") ;
            hold on
            axis equal
            plot([min(dhdB) max(dhdB)],[min(dhdB) max(dhdB)],"--k") ;
            axis equal tight ;
            xlabel(" $d\dot{h}/dB$",Interpreter="latex")  ;
            ylabel("Finite difference $d\dot{h}/dB$",Interpreter="latex")
            ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
            axis on ; axis equal tight ; box off
            title("Comparision between adjoint and finite-differences gradient calculations")
            set(gcf,'Color','white')

        end



    end

    %% C
    if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
        % du/dC
        dudC=dudC(:,NodeTest);
        dvdC=dvdC(:,NodeTest);
        if ~isempty(dhdC)
            dhdC=dhdC(:,NodeTest);
        end
        %%
        F=Funperturbed;

        C0=F.C;
        dC=F.C(NodeTest)*1e-1;

        Cp=C0;
        Cp(NodeTest)=Cp(NodeTest)+dC;

        F.C=Cp;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
        up=F.ub; vp=F.vb; dhdtp=dhdt;

        F.C=C0;
        Cm=C0;
        Cm(NodeTest)=Cm(NodeTest)-dC;
        F.C=Cm;
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
        um=F.ub; vm=F.vb; dhdtm=dhdt;

        F.C=C0;

        dudCpert=(up-um)/(2*dC) ;
        dvdCpert=(vp-vm)/(2*dC) ;
        dhdotdCpert=(dhdtp-dhdtm)/(2*dC) ;

        scale=log(10)*C0;
        dudCpert=dudCpert.*scale ;
        dvdCpert=dvdCpert.*scale ;
        dhdotdCpert=dhdotdCpert.*scale ;


        % dv/dC
        figCu=FindOrCreateFigure("du/dC comparision"); clf(figCu)
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudC,CreateNewFigure=false)  ; title("$du/dC$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudCpert,CreateNewFigure=false) ; title("$du/dC$ finite differences ($\log$ scale)",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,dudC-dudCpert,CreateNewFigure=false) ; title("$du/dC$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";

        figCv=FindOrCreateFigure("dv/dC comparision"); clf(figCv)
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdC,CreateNewFigure=false)  ; title("$dv/dC$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdCpert,CreateNewFigure=false) ; title("$dv/dA$ finite differences ($\log$ scale)",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,dvdC-dvdCpert,CreateNewFigure=false) ; title("$dv/dA$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";

        if ~isempty(dhdC)


            figChdot=FindOrCreateFigure("dhdot/dC comparision"); clf(figChdot)
            T=tiledlayout("flow");

            T1=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdC,CreateNewFigure=false)  ; title("$d\dot{h}/dC$ sensitvity ",Interpreter="latex") ; subtitle("")
            title(cbar,"")

            T2=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdotdCpert,CreateNewFigure=false) ; title("$d\dot{h}/dC$ finite differences ",Interpreter="latex") ; subtitle("") ; title(cbar,"")

            T3=nexttile;
            UaPlots(CtrlVar,MUA,F,dhdC-dhdotdCpert,CreateNewFigure=false) ; title("$d\dot{h}/dC$ differences",Interpreter="latex") ; subtitle("")
            CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

            T.Padding="loose";   T.TileSpacing="tight";


        end

        figCgradu=FindOrCreateFigure("du/dC gradient test") ;  clf(figCgradu)
        plot(dudC,dudCpert,"or") ;
        hold on
        axis equal
        AX=axis;
        plot([min(dudC) max(dudC)],[min(dudC) max(dudC)],"--k") ;
        axis equal tight ;
        xlabel(" $du/dC$",Interpreter="latex")  ;
        ylabel("Finite difference $du/dC$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')

        figCgradv=FindOrCreateFigure("dv/dC gradient test") ;  clf(figCgradv)
        plot(dvdC,dvdCpert,"or") ;
        hold on
        axis equal
        plot([min(dvdC) max(dvdC)],[min(dvdC) max(dvdC)],"--k") ;
        axis equal tight ;
        xlabel(" $dv/dC$",Interpreter="latex")  ;
        ylabel("Finite difference $dv/dC$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')



        if ~isempty(dhdC)

            figCgradh=FindOrCreateFigure("dhdot/dC gradient test") ;  clf(figCgradh)
            plot(dhdC,dhdotdCpert,"or") ;
            hold on
            axis equal
            plot([min(dhdC) max(dhdC)],[min(dhdC) max(dhdC)],"--k") ;
            axis equal tight ;
            xlabel(" $d\dot{h}/dC$",Interpreter="latex")  ;
            ylabel("Finite difference $d\dot{h}/dC$",Interpreter="latex")
            ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
            axis on ; axis equal tight ; box off
            title("Comparision between adjoint and finite-differences gradient calculations")
            set(gcf,'Color','white')

        end


        %%
    end


end