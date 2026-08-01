
function [dudA,dvdA,dhdA,dudB,dvdB,dhdB,dudC,dvdC,dhdC]=duv_hdABC(UserVar,CtrlVar,RunInfo,MUA,F,l,BCs)


narginchk(7,7)
nargoutchk(9,9)

dudA=[]; dvdA=[]; dhdA=[];
dudB=[]; dvdB=[]; dhdB=[];
dudC=[]; dvdC=[]; dhdC=[];

[UserVar,RunInfo,F,l,KdFuvduv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);

% CtrlVar.uvAssembly.ZeroFields=false;
% CtrlVar.uvMatrixAssembly.Ronly=false;
% [~,KdFuvduv]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs);


if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)
    tA=tic;
    [dudA,dvdA,dhdA]=duvhdotdAFunc(CtrlVar,MUA,F,l,BCs,KdFuvduv) ;  % this has been tested against finite-differences and is good, also for dhdotdA
    tA=toc(tA);
    fprintf("A sensitivities for %i nodes calculated in %f sec\n",MUA.Nnodes,tA)
    scale=log(10)*F.AGlen;
    dudA=dudA.*scale ; % using implicit expansion
    dvdA=dvdA.*scale ; % using implicit expansion
    if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
        dhdA=dhdA.*scale ;
    else
        dhdA=[];
    end
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

    scale=log(10)*F.C;  % this has to be a row vector
    dudC=dudC.*scale ;   % using implicit expansion
    dvdC=dvdC.*scale ;   % using implicit expansion
    if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
        dhdC=dhdC.*scale ;
    else
        dhdC=[];
    end
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




TestSensitivites=true;

if TestSensitivites

    %% Test


    Funperturbed=F;
    NodeTest=804;
    NodeTest=1200;
    %NodeTest=1500;

    %% A
    if contains(CtrlVar.Inverse.InvertFor,"logaglen",IgnoreCase=true)


        NodeTest=2500;
        F=Funperturbed;
        DeltaRel=0.01;
        Field="AGlen";
        [dudApert,dvdApert,dhdotdApert]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

        PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudA,dvdA,dhdA,dudApert,dvdApert,dhdotdApert);


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

        NodeTest=2500;
        F=Funperturbed;
        DeltaRel=0.1;
        Field="C";
        [dudCpert,dvdCpert,dhdotdCpert]=FiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel);

        PlotModelAndFiniteDifferenceSensitivities(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Field,NodeTest,DeltaRel,dudC,dvdC,dhdC,dudCpert,dvdCpert,dhdotdCpert);





        %%
    end


end