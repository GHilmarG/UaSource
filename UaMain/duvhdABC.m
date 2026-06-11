
% function [KdudA,KdvdA,KdudB,KdvdB,KdudC,KdvdC]=duvhdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs)

function   [KdudA,KdvdA,KdhdA,Kdudh,Kdvdh,Kdhdh,KdudC,KdvdC,KdhdC,KdFuvhduvh]=duvhdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs)

narginchk(7,7)

KdudA=[]; KdvdA=[];
Kdudh=[]; Kdvdh=[];
KdudC=[]; KdvdC=[];

KdFuvhduvh=[];


if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    tA=tic;
    [KdudA,KdvdA,KdhdA,KdFuvhduvh]=duvhdA(CtrlVar,MUA,F,l,BCs,KdFuvhduvh) ;

    tA=toc(tA);
    fprintf("dudA and dvdA sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    tB=tic;

    % [KdudB,KdvdB]=duvdBFunc(CtrlVar,MUA,F,l,BCs) ;  % this has been tested against finite-differences and is good
    [Kdudh,Kdvdh,Kdhdh,KdFuvhduvh]=duvhdh(CtrlVar,MUA,F,l,BCs,KdFuvhduvh) ;
    tB=toc(tB);
    fprintf("dudh and dvdh sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tB)
end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
    tC=tic;

    %[KdudC,KdvdC]=duvdCFunc(CtrlVar,MUA,F,l,BCs) ; % this has been tested against finite-differences and is good
    [KdudC,KdvdC,KdhdC,KdFuvhduvh]=duvhdC(CtrlVar,MUA,F,l,BCs,KdFuvhduvh) ;

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
    KdudA=KdudA.*scale ; % using implicit expansion
    KdvdA=KdvdA.*scale ; % using implicit expansion
    KdhdA=KdhdA.*scale ;

end

if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)

    scale=log(10)*F.C';  % this has to be a row vector
    KdudC=KdudC.*scale ;   % using implicit expansion
    KdvdC=KdvdC.*scale ;   % using implicit expansion
    KdhdC=KdhdC.*scale ;

end

TestSensitivites=true;

if TestSensitivites

    %% Test



    NodeTest=804;
    %NodeTest=1200;
    %NodeTest=1500;

    %% A
    if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)


        % du/dA  : these are now log10 sensitivities
        KdudA=KdudA(:,NodeTest);
        KdvdA=KdvdA(:,NodeTest);

        if ~isempty(KdhdA)
            KdhdA=KdhdA(:,NodeTest);
        end

        %%

        % solve the diagnostic problem
        [UserVar,RunInfo,F0,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        F1=F0;

        A0=F0.AGlen;
        dA=F0.AGlen(NodeTest)*0.0001;

        Ap=A0;
        Ap(NodeTest)=Ap(NodeTest)+dA;
        F1.AGlen=Ap;
        [UserVar,RunInfo,F1,l,BCs]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l,l,BCs);
        up=F1.ub; vp=F1.vb; dhdtp=F1.dhdt;

        F1=F0;
        Am=A0;
        Am(NodeTest)=Am(NodeTest)-dA;
        F1.AGlen=Am;
        [UserVar,RunInfo,F1,l,BCs]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l,l,BCs);
        um=F1.ub; vm=F1.vb; dhdtm=F1.dhdt;



        dudApert=(up-um)/(2*dA) ;
        dvdApert=(vp-vm)/(2*dA) ;
        dhdApert=(dhdtp-dhdtm)/(2*dA) ;

        scale=log(10)*A0;
        dudApert=dudApert.*scale ;
        dvdApert=dvdApert.*scale ;
        dhdApert=dhdApert.*scale ;

        % dv/dA
        figAu=FindOrCreateFigure("du/dA comparision");
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,KdudA,CreateNewFigure=false)  ; title("$du/dA$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudApert,CreateNewFigure=false) ; title("$du/dA$ finite differences ($\log$ scale) ",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,KdudA-dudApert,CreateNewFigure=false) ; title("$du/dA$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";

        figAv=FindOrCreateFigure("dv/dA comparision"); clf(figAv)
        T=tiledlayout("flow");

        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,KdvdA,CreateNewFigure=false)  ; title("$dv/dA$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")

        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdApert,CreateNewFigure=false) ; title("$dv/dA$ finite differences ($\log$ scale)",Interpreter="latex") ; subtitle("") ; title(cbar,"")

        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,KdvdA-dvdApert,CreateNewFigure=false) ; title("$dv/dA$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

        T.Padding="loose";   T.TileSpacing="tight";

        % dh sensitivity to A

        if ~isempty(KdhdA)


            % dh/dA
            figAh=FindOrCreateFigure("dh/dA comparision"); clf(figAh)
            T=tiledlayout("flow");

            T1=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,KdhdA,CreateNewFigure=false)  ; title("$dh/dA$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
            title(cbar,"")

            T2=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdApert,CreateNewFigure=false) ; title("$dh/dA$ finite differences ($\log$ scale) ",Interpreter="latex") ; subtitle("") ; title(cbar,"")

            T3=nexttile;
            UaPlots(CtrlVar,MUA,F,KdhdA-dhdApert,CreateNewFigure=false) ; title("$dh/dA$ differences",Interpreter="latex") ; subtitle("")
            CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

            T.Padding="loose";   T.TileSpacing="tight";


        end



        figAgradu=FindOrCreateFigure("du/dA gradient test") ;  clf(figAgradu)
        plot(KdudA,dudApert,"or") ;
        hold on
        axis equal
        AX=axis;
        plot([min(KdudA) max(KdudA)],[min(KdudA) max(KdudA)],"--k") ;
        axis equal tight ;
        xlabel(" $du/dA$",Interpreter="latex")  ;
        ylabel("Finite difference $du/dA$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')

        figAgradv=FindOrCreateFigure("dv/dA gradient test") ;  clf(figAgradv)
        plot(KdvdA,dvdApert,"or") ;
        hold on
        axis equal
        plot([min(KdvdA) max(KdvdA)],[min(KdvdA) max(KdvdA)],"--k") ;
        axis equal tight ;
        xlabel(" $dv/dA$",Interpreter="latex")  ;
        ylabel("Finite difference $dv/dA$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')


        if ~isempty(KdhdA)

            figAgradh=FindOrCreateFigure("dh/dA gradient test") ;  clf(figAgradh)
            plot(KdhdA,dhdApert,"or") ;
            hold on
            axis equal
            plot([min(KdhdA) max(KdhdA)],[min(KdhdA) max(KdhdA)],"--k") ;
            axis equal tight ;
            xlabel(" $dh/dA$",Interpreter="latex")  ;
            ylabel("Finite difference $dh/dA$",Interpreter="latex")
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

        CtrlVar.Calculate.Geometry="bh-FROM-sBS" ;

        Kdudh=Kdudh(:,NodeTest);
        Kdvdh=Kdvdh(:,NodeTest);

        if ~isempty(Kdhdh)
            Kdhdh=Kdhdh(:,NodeTest);
        end

        % solve the unperturbed diagnostic problem
        [UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
    
        
        F0=F; 
        dh=F0.h(NodeTest)*0.01;

        hp=F0.h;
        hp(NodeTest)=hp(NodeTest)+dh;
        F0.h=hp;
        [F0.b,F0.s,F0.h,F0.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F0.h,F0.S,F0.B,F0.rho,F0.rhow,F0.GF);

        [UserVar,RunInfo,F0,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,l);
        [UserVar,RunInfo,F1,l,BCs]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F0,l,l,BCs);
        up=F0.ub; vp=F0.vb; dhdtp=(F1.h-F0.h)/CtrlVar.dt;


        F0=F;
        hm=F0.h;
        hm(NodeTest)=hm(NodeTest)-dh;
        F0.h=hm;
        [UserVar,RunInfo,F0,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,l);
        [UserVar,RunInfo,F1,l,BCs]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F0,l,l,BCs);
        um=F0.ub; vm=F0.vb; dhdtm=(F1.h-F0.h)/CtrlVar.dt;



        dudhpert=(up-um)/(2*dh) ;
        dvdhpert=(vp-vm)/(2*dh) ;
        dhdhpert=(dhdtp-dhdtm)/(2*dh) ;



        % dv/dB
        figBu=FindOrCreateFigure("du/dh comparision");
        T=tiledlayout("flow");
        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,Kdudh,CreateNewFigure=false)  ; title("$du/dh$ sensitvity",Interpreter="latex") ; subtitle("")
        title(cbar,"")
        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudhpert,CreateNewFigure=false) ; title("$du/dh$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")
        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,Kdudh-dudhpert,CreateNewFigure=false) ; title("$du/dh$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
        T.Padding="loose";   T.TileSpacing="tight";



        figBv=FindOrCreateFigure("dv/dh comparision");
        T=tiledlayout("flow");
        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,Kdvdh,CreateNewFigure=false)  ; title("$dv/dh$ sensitvity",Interpreter="latex") ; subtitle("")
        title(cbar,"")
        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdhpert,CreateNewFigure=false) ; title("$dv/dh$ finite differences",Interpreter="latex") ; subtitle("") ; title(cbar,"")
        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,Kdvdh-dvdhpert,CreateNewFigure=false) ; title("$dv/dh$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
        T.Padding="loose";   T.TileSpacing="tight";



        if ~isempty(Kdhdh)

            % dh/dh
            fighh=FindOrCreateFigure("dh/dh comparision"); clf(fighh)
            T=tiledlayout("flow");
            T1=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,Kdhdh,CreateNewFigure=false)  ; title("$dh/dh$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
            title(cbar,"")
            T2=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdhpert,CreateNewFigure=false) ; title("$dh/dh$ finite differences ($\log$ scale) ",Interpreter="latex") ; subtitle("") ; title(cbar,"")
            T3=nexttile;
            UaPlots(CtrlVar,MUA,F,Kdhdh-dhdhpert,CreateNewFigure=false) ; title("$dh/dh$ differences",Interpreter="latex") ; subtitle("")
            CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
            T.Padding="loose";   T.TileSpacing="tight";

        end



        figBgrad=FindOrCreateFigure("du/dh gradient test") ;  clf(figBgrad)
        plot(Kdudh,dudhpert,"or") ;
        hold on
        axis equal
        AX=axis;
        plot([min(Kdudh) max(Kdudh)],[min(Kdudh) max(Kdudh)],"--k") ;
        axis equal tight ;
        xlabel(" $du/dh$",Interpreter="latex")  ;
        ylabel("Finite difference $du/dh$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')

        figBgrad=FindOrCreateFigure("dv/dh gradient test") ;  clf(figBgrad)
        plot(Kdvdh,dvdhpert,"or") ;
        hold on
        axis equal
        plot([min(Kdvdh) max(Kdvdh)],[min(Kdvdh) max(Kdvdh)],"--k") ;
        axis equal tight ;
        xlabel(" $dv/dh$",Interpreter="latex")  ;
        ylabel("Finite difference $dv/dh$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')


        if ~isempty(Kdhdh)

            fighgradh=FindOrCreateFigure("dh/dh gradient test") ;  clf(fighgradh)
            plot(Kdhdh,dhdhpert,"or") ;
            hold on
            axis equal
            plot([min(Kdhdh) max(Kdhdh)],[min(Kdhdh) max(Kdhdh)],"--k") ;
            axis equal tight ;
            xlabel(" $dh/dh$",Interpreter="latex")  ;
            ylabel("Finite difference $dh/dh$",Interpreter="latex")
            ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
            axis on ; axis equal tight ; box off
            title("Comparision between adjoint and finite-differences gradient calculations")
            set(gcf,'Color','white')

        end



    end

    %% C
    if contains(CtrlVar.Inverse.InvertFor,"logc",IgnoreCase=true)
        % du/dC
        KdudC=KdudC(:,NodeTest);
        KdvdC=KdvdC(:,NodeTest);
        if ~isempty(KdhdC)
            KdhdC=KdhdC(:,NodeTest);
        end
        %%

        % solve the diagnostic problem
        [UserVar,RunInfo,F0,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        F1=F0;

        C0=F0.C;
        dC=F0.C(NodeTest)*0.0001;

        Cp=C0;
        Cp(NodeTest)=Cp(NodeTest)+dC;
        F1.C=Cp;
        [UserVar,RunInfo,F1,l,BCs]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l,l,BCs);
        up=F1.ub; vp=F1.vb; dhdtp=F1.dhdt;

        F1=F0;
        Cm=C0;
        Cm(NodeTest)=Cm(NodeTest)-dC;
        F1.C=Cm;
        [UserVar,RunInfo,F1,l,BCs]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l,l,BCs);
        um=F1.ub; vm=F1.vb; dhdtm=F1.dhdt;


        dudCpert=(up-um)/(2*dC) ;
        dvdCpert=(vp-vm)/(2*dC) ;
        dhdCpert=(dhdtp-dhdtm)/(2*dC) ;

        scale=log(10)*C0;
        dudCpert=dudCpert.*scale ;
        dvdCpert=dvdCpert.*scale ;
        dhdCpert=dhdCpert.*scale ;


        % dv/dC
        figCu=FindOrCreateFigure("du/dC comparision"); clf(figCu)
        T=tiledlayout("flow");
        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,KdudC,CreateNewFigure=false)  ; title("$du/dC$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")
        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dudCpert,CreateNewFigure=false) ; title("$du/dC$ finite differences ($\log$ scale)",Interpreter="latex") ; subtitle("") ; title(cbar,"")
        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,KdudC-dudCpert,CreateNewFigure=false) ; title("$du/dC$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
        T.Padding="loose";   T.TileSpacing="tight";

        figCv=FindOrCreateFigure("dv/dC comparision"); clf(figCv)
        T=tiledlayout("flow");
        T1=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,KdvdC,CreateNewFigure=false)  ; title("$dv/dC$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
        title(cbar,"")
        T2=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,dvdCpert,CreateNewFigure=false) ; title("$dv/dA$ finite differences ($\log$ scale)",Interpreter="latex") ; subtitle("") ; title(cbar,"")
        T3=nexttile;
        UaPlots(CtrlVar,MUA,F,KdvdC-dvdCpert,CreateNewFigure=false) ; title("$dv/dA$ differences",Interpreter="latex") ; subtitle("")
        CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
        T.Padding="loose";   T.TileSpacing="tight";


        % dh sensitivity to C

        if ~isempty(KdhdC)


            % dh/dC
            figCh=FindOrCreateFigure("dh/dC comparision"); clf(figCh)
            T=tiledlayout("flow");
            T1=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,KdhdC,CreateNewFigure=false)  ; title("$dh/dC$ sensitvity ($\log$ scale)",Interpreter="latex") ; subtitle("")
            title(cbar,"")
            T2=nexttile;
            cbar=UaPlots(CtrlVar,MUA,F,dhdCpert,CreateNewFigure=false) ; title("$dh/dC$ finite differences ($\log$ scale) ",Interpreter="latex") ; subtitle("") ; title(cbar,"")
            T3=nexttile;
            UaPlots(CtrlVar,MUA,F,KdhdC-dhdCpert,CreateNewFigure=false) ; title("$dh/dC$ differences",Interpreter="latex") ; subtitle("")
            CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
            T.Padding="loose";   T.TileSpacing="tight";

        end


        figCgradu=FindOrCreateFigure("du/dC gradient test") ;  clf(figCgradu)
        plot(KdudC,dudCpert,"or") ;
        hold on
        axis equal
        AX=axis;
        plot([min(KdudC) max(KdudC)],[min(KdudC) max(KdudC)],"--k") ;
        axis equal tight ;
        xlabel(" $du/dC$",Interpreter="latex")  ;
        ylabel("Finite difference $du/dC$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')

        figCgradv=FindOrCreateFigure("dv/dC gradient test") ;  clf(figCgradv)
        plot(KdvdC,dvdCpert,"or") ;
        hold on
        axis equal
        plot([min(KdvdC) max(KdvdC)],[min(KdvdC) max(KdvdC)],"--k") ;
        axis equal tight ;
        xlabel(" $dv/dC$",Interpreter="latex")  ;
        ylabel("Finite difference $dv/dC$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')


        if ~isempty(KdhdC)

            figCgradh=FindOrCreateFigure("dh/dC gradient test") ;  clf(figCgradh)
            plot(KdhdC,dhdCpert,"or") ;
            hold on
            axis equal
            plot([min(KdhdC) max(KdhdC)],[min(KdhdC) max(KdhdC)],"--k") ;
            axis equal tight ;
            xlabel(" $dh/dC$",Interpreter="latex")  ;
            ylabel("Finite difference $dh/dC$",Interpreter="latex")
            ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
            axis on ; axis equal tight ; box off
            title("Comparision between adjoint and finite-differences gradient calculations")
            set(gcf,'Color','white')

        end



        %%
    end

end