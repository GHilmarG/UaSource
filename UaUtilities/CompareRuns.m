function [Figs,F]=CompareRuns(files,Var)


persistent FigCounter


if isempty(FigCounter) ; FigCounter=1; end

N=numel(files);

if ischar(files{1})
    
    
    for I=1:N
        F{I}=load(files{I});
    end
    
    for I=1:N
        if ~isfield(F{I},'CtrlVar')  && isfield(F{I},'CtrlVarInRestartFile')
            F{I}.CtrlVar=F{I}.CtrlVarInRestartFile;
        end
    end
    
    %     if N==2
    %         F{3}=F{1} ;
    %         F{3}.s=F{1}.s-F{2}.s;
    %         F{3}.b=F{1}.b-F{2}.b;
    %         F{3}.h=F{1}.h-F{2}.h;
    %         F{3}.B=F{1}.B-F{2}.B;
    %         F{3}.ub=F{1}.ub-F{2}.ub;
    %         F{3}.vb=F{1}.vb-F{2}.vb;
    %         F{3}.dhdt=F{1}.dhdt-F{2}.dhdt;
    %         N=N+1;
    %     end
    
    
    for I=1:N
        if ~isfield(F{I},'MUA')
            F{I}.CtrlVar.niph=F{I}.nip ;
            F{I}.CtrlVar.nip=F{I}.nip ;
            F{I}.MUA=CreateMUA(F{I}.CtrlVar,F{I}.connectivity,F{I}.coordinates);
        end
        
        if ~isfield(F{I},'GF')
            F{I}.GF = GL2d(F{I}.B,F{I}.S,F{I}.h,F{I}.rhow,F{I}.rho,F{I}.MUA.connectivity,F{I}.CtrlVar);
        end
        
        if ~isfield(F{I},'ub')
            F{I}.ub=F{I}.u;
            F{I}.vb=F{I}.v;
        end
        
    end
    
    
    for I=1:N
        F{I}.GLgeo=GLgeometry(F{I}.MUA.connectivity,F{I}.MUA.coordinates,F{I}.GF,F{I}.CtrlVar);
    end
else
    F=files;
end

if strcmp(Var,'ubvb')
    for I=1:N
        Figs{I}=figure(FigCounter); FigCounter=FigCounter+1;
        k=2;
        
        PlotBoundary(F{I}.MUA.Boundary,F{I}.MUA.connectivity,F{I}.MUA.coordinates,F{I}.CtrlVar,'b')
        hold on
        F{I}.CtrlVar.VelPlotIntervalSpacing='log10';
        QuiverColorGHG(F{I}.MUA.coordinates(1:k:end,1),F{I}.MUA.coordinates(1:k:end,2),F{I}.ub(1:k:end),F{I}.vb(1:k:end),F{I}.CtrlVar);
        if I==3
            hold on
            plot(F{1}.GLgeo(:,[3 4])'/F{1}.CtrlVar.PlotXYscale,F{1}.GLgeo(:,[5 6])'/F{1}.CtrlVar.PlotXYscale,'k','LineWidth',1);
            plot(F{2}.GLgeo(:,[3 4])'/F{2}.CtrlVar.PlotXYscale,F{2}.GLgeo(:,[5 6])'/F{2}.CtrlVar.PlotXYscale,'r','LineWidth',1);
        else
            hold on ;
            plot(F{I}.GLgeo(:,[3 4])'/F{I}.CtrlVar.PlotXYscale,F{I}.GLgeo(:,[5 6])'/F{1}.CtrlVar.PlotXYscale,'k','LineWidth',1);
        end
    end
elseif strcmp(Var,'mesh')
    for I=1:N
        
        Figs{I}=figure(FigCounter); FigCounter=FigCounter+1;
        PlotFEmesh(F{I}.MUA.coordinates,F{I}.MUA.connectivity,F{I}.CtrlVar)
        title(sprintf('#Ele=%-i, #Nodes=%-i, #nod=%-i',...
            F{I}.MUA.Nele,F{I}.MUA.Nnodes,F{I}.MUA.nod))
        hold on ; plot(F{I}.GLgeo(:,[3 4])'/F{I}.CtrlVar.PlotXYscale,F{I}.GLgeo(:,[5 6])'/F{I}.CtrlVar.PlotXYscale,'r','LineWidth',1);
    end
else
    
    for I=1:N
        
        Figs{I}=figure(FigCounter); FigCounter=FigCounter+1;
        PlotMeshScalarVariable(F{I}.CtrlVar,F{I}.MUA,eval(['F{',num2str(I),'}.',Var]));
        F{I}.GLgeo=GLgeometry(F{I}.MUA.connectivity,F{I}.MUA.coordinates,F{I}.GF,F{I}.CtrlVar);
        if I==3
            hold on
            plot(F{1}.GLgeo(:,[3 4])'/F{1}.CtrlVar.PlotXYscale,F{1}.GLgeo(:,[5 6])'/F{1}.CtrlVar.PlotXYscale,'k','LineWidth',1);
            plot(F{2}.GLgeo(:,[3 4])'/F{2}.CtrlVar.PlotXYscale,F{2}.GLgeo(:,[5 6])'/F{2}.CtrlVar.PlotXYscale,'r','LineWidth',1);
        else
            hold on ;
            plot(F{I}.GLgeo(:,[3 4])'/F{I}.CtrlVar.PlotXYscale,F{I}.GLgeo(:,[5 6])'/F{1}.CtrlVar.PlotXYscale,'k','LineWidth',1);
        end
        
    end
    
end


end
