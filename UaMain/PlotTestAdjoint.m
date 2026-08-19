
function PlotTestAdjoint(CtrlVar,MUA,F,InvFinalValues)

%


if ~isempty(InvFinalValues.dJdAGlenTest)
    IA=find(~isnan(InvFinalValues.dJdAGlenTest)) ;
    fprintf('------------------------------------ dJ/dA  ---------------------------------------------------------------------\n')
    fprintf('#Node/Ele  dJdA          dJdATest      dJdA-dJdATest     dJdA/dtdATest  (dJdA-dJdATest)/dJdA \n')

    for ii=1:numel(IA)
        I=IA(ii);
        fprintf('%i %15g %15g  %15g  %15g %15g \n',I,...
            InvFinalValues.dJdAGlen(I),...
            InvFinalValues.dJdAGlenTest(I),...
            InvFinalValues.dJdAGlen(I)-InvFinalValues.dJdAGlenTest(I),...
            InvFinalValues.dJdAGlen(I)/InvFinalValues.dJdAGlenTest(I),...
            (InvFinalValues.dJdAGlen(I)-InvFinalValues.dJdAGlenTest(I))/InvFinalValues.dJdAGlen(I))
    end

    figAgrad=FindOrCreateFigure("dJ/dA test") ;  clf(figAgrad)
    plot(InvFinalValues.dJdAGlen,InvFinalValues.dJdAGlenTest,"or") ;
    hold on
    plot(InvFinalValues.dJdAGlen,InvFinalValues.dJdAGlen,"--k") ;

    xlabel("Adjoint $dJ/dA$",Interpreter="latex")  ;
    ylabel("Finite difference $dJ/dA$",Interpreter="latex")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    axis on ; axis equal tight ; box off
    %       axis([min(InvFinalValues.dJACTest) max(InvFinalValues.dJACTest) min(InvFinalValues.dJdATest) max(InvFinalValues.dJdATest)])
    title("Comparision betweenadjoint and finite-differences gradient calculations")
    set(gcf,'Color','white')



end

if ~isempty(InvFinalValues.dJdCTest)

    IC=find(~isnan(InvFinalValues.dJdCTest)) ;

    fprintf('--------------------------------------- dJ/dC ----------------------------------------------------------------------\n')

    fprintf('#Node/Ele  dJdC          dJdCTest      dJdC-dJdCTest     dJdC/dtdCTest   (dJdC-dJdCTest)/dJdC\n')

    for ii=1:numel(IC)
        I=IC(ii);
        fprintf('%i %15g %15g  %15g  %15g %15g \n',I,...
            InvFinalValues.dJdC(I),...
            InvFinalValues.dJdCTest(I),...
            InvFinalValues.dJdC(I)-InvFinalValues.dJdCTest(I),...
            InvFinalValues.dJdC(I)/InvFinalValues.dJdCTest(I),...
            (InvFinalValues.dJdC(I)-InvFinalValues.dJdCTest(I))/InvFinalValues.dJdC(I))
    end

    %%

    figCgrad=FindOrCreateFigure("dJ/dC test") ;  clf(figCgrad)
    plot(InvFinalValues.dJdC,InvFinalValues.dJdCTest,"or") ;
    hold on
    plot(InvFinalValues.dJdC,InvFinalValues.dJdC,"--k") ;
    xlabel("Adjoint $dJ/dC$",Interpreter="latex")  ;
    ylabel("Finite difference $dJ/dC$",Interpreter="latex")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    axis on ; axis equal tight;
    %       axis([min(InvFinalValues.dJdCTest) max(InvFinalValues.dJdCTest) min(InvFinalValues.dJdCTest) max(InvFinalValues.dJdCTest)])
    box off
    title("Comparision betweenadjoint and finite-differences gradient calculations")
    set(gcf,'Color','white')
end

if ~isempty(InvFinalValues.dJdBTest)
    fprintf('--------------------------------------- dJ/dB ----------------------------------------------------------------------\n')

    IB=find(~isnan(InvFinalValues.dJdBTest)) ;

    fprintf('#Node/Ele  dJdB          dJdBTest      dJdB-dJdBTest     dJdB/dtdBTest   \n')

    for ii=1:numel(IB)
        I=IB(ii);
        fprintf('%i %15g %15g  %15g  %15g %15g \n',I,...
            InvFinalValues.dJdB(I),...
            InvFinalValues.dJdBTest(I),...
            InvFinalValues.dJdB(I)-InvFinalValues.dJdBTest(I),...
            InvFinalValues.dJdB(I)/InvFinalValues.dJdBTest(I),...
            (InvFinalValues.dJdB(I)-InvFinalValues.dJdBTest(I))/InvFinalValues.dJdB(I))
    end


    fprintf('--------------------------------------------------------------------------------------------------------------------------\n')

    %%

    figBgrad=FindOrCreateFigure("dJ/dB test") ;  clf(figBgrad)
    plot(InvFinalValues.dJdB,InvFinalValues.dJdBTest,"or") ;
    hold on
    plot(InvFinalValues.dJdB,InvFinalValues.dJdB,"--k") ;
    axis equal tight ;
    axis([min(InvFinalValues.dJdBTest) max(InvFinalValues.dJdBTest) min(InvFinalValues.dJdBTest) max(InvFinalValues.dJdBTest)])
    box off
    xlabel("Adjoint $dJ/dB$",Interpreter="latex")  ;
    ylabel("Finite difference $dJ/dB$",Interpreter="latex")
    ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
    title("Comparision between adjoint and finite-differences gradient calculations")
    set(gcf,'Color','white')
end
%%
%[dJdp(iRange) dJdpTest(iRange)   dJdp(iRange)-dJdpTest(iRange) dJdp(iRange)./dJdpTest(iRange)]
% iRange=find(~isnan(dJdpTest));
% dJdpTest(iRange)-dJdp(iRange))/norm(dJdp(iRange)))
% fprintf('Norm test: ||dJdpTest-dJdp||/||dJdp||= %g \n ',norm(dJdpTest(iRange)-dJdp(iRange))/norm(dJdp(iRange)))

%%

if ~(isempty(InvFinalValues.dJdC) && isempty(InvFinalValues.dJdCTest))

    IFigC=FindOrCreateFigure("dJ/dC test over mesh") ; clf(IFigC);

    TileC=tiledlayout("flow");

    TdJdC=nexttile;

    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('dJdC Adjoint directional derivative')

    TdJdCTest=nexttile;
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdCTest) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('$dJ/dC$ Brute force directional derivative','interpreter','latex')


    TdJdCDiff=nexttile;
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC-InvFinalValues.dJdCTest) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('Difference between adjoint and brute force directional derivatives')


    TdJdCRatio=nexttile;
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC./InvFinalValues.dJdCTest) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('Ratio between adjoint and brute force directional derivatives')

    TileC.TileSpacing='tight';
    TileC.Padding='tight';
    %IFigC.Position=[948.43 41.571 1246.3 1115.4];
    %%
end


if ~(isempty(InvFinalValues.dJdAGlen) && isempty(InvFinalValues.dJdAGlenTest))


    IFigAGlen=FindOrCreateFigure("dJ/dA test over mesh") ; clf(IFigAGlen);

    TileA=tiledlayout("flow");

    nexttile
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('dJdAGlen Adjoint directional derivative')

    nexttile
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlenTest) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('dJdAGlen Brute force directional derivative')


    nexttile
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen-InvFinalValues.dJdAGlenTest) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('Difference between adjoint and brute force directional derivatives')

    nexttile
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen./InvFinalValues.dJdAGlenTest) ;
    hold on
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    title('Ratio between adjoint and brute force directional derivatives')

    TileA.TileSpacing='tight';
    TileA.Padding='tight';
    % IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
    %%
end


if ~(isempty(InvFinalValues.dJdB) && isempty(InvFinalValues.dJdBTest))

    %%

    IFigb=FindOrCreateFigure("dJ/dB test over mesh") ; clf(IFigb)
    TileB=tiledlayout("flow");
    TdJdB=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB,PlotUnderMesh=true,CreateNewFigure=false);
    title("$dJ/dB$ Adjoint directional derivative")
    subtitle("")

    TdJdBTest=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false);
    title('$dJ/dB$ Brute force directional derivative',Interpreter='latex')

    TdJdBDiff=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB-InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false);
    title('Difference between adjoint and brute force derivatives')
    subtitle("")
    axis(TdJdBDiff) ; CM=cmocean('balanced',25,'pivot',0) ; colormap(TdJdBDiff,CM);

    TdJdBRatio=nexttile;
    UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB./InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false) ;
    title('Ratio between adjoint and brute force derivatives')
    subtitle("")

    %  IFigb.Position=[1.5714 41.571 1096 1115.4];
    TileB.TileSpacing='tight';
    TileB.Padding='tight';
    %CL=TdJdB.CLim; TdJdBTest.CLim=CL;  TdJdBDiff.CLim=CL; TdJdBRatio.CLim=[0.7 1.3];
    axis(TdJdBDiff) ; CM=cmocean('balanced',25,'pivot',0) ; colormap(TdJdBDiff,CM);
    %%
end

