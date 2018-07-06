function UserVar=UaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%%
% This routine is called during the run and can be used for saving and/or plotting data.
% 
%   UserVar=UaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%
% Write your own version of this routine and put it in you local run directory.
% 
% Inputs:
% 
%   BCs             Structure with all boundary conditions
%   UaVars          A structure with fields such as
%                   s, b, S, B, rho, rhow, ub, vb, ud, vd, AGlen, n , C, m , as, ab, g
%   l               Lagrange parameters related to the enforcement of boundary
%                   conditions.
%   GF              Grounding floating mask for nodes and elements.
%
% If preferred to work directly with the variables rather than the respective fields of the structure 
% `UaVars', then`UaVars' can easily be converted into variables using v2struc.
%
%%
v2struct(F);
time=CtrlVar.time; 
plots='-plot-';

if contains(plots,'-save-')
    
    % save data in files with running names
    % check if folder 'ResultsFiles' exists, if not create
    
    if exist(fullfile(cd,UserVar.Outputsdirectory),'dir')~=7
        mkdir(CtrlVar.Outputsdirectory) ;
    end
    
    if strcmp(CtrlVar.UaOutputsInfostring,'Last call')==0
  
        
        FileName=sprintf('%s/%07i-Nodes%i-Ele%i-Tri%i-kH%i-%s.mat',...
            CtrlVar.Outputsdirectory,round(100*time),MUA.Nnodes,MUA.Nele,MUA.nod,1000*CtrlVar.kH,CtrlVar.Experiment);
        fprintf(' Saving data in %s \n',FileName)
        save(FileName,'UserVar','CtrlVar','MUA','time','s','b','S','B','h','ub','vb','C','dhdt','AGlen','m','n','rho','rhow','as','ab','GF')
        
    end
    
end

if contains(plots,'-plot-')
    
    figsWidth=1000 ; figHeights=300;
    GLgeo=[]; xGL=[] ; yGL=[];
    %%
    fig100=figure(100) ;
    fig100.Position=[50 50 figsWidth 3*figHeights];
    subplot(4,1,1)
    PlotMeshScalarVariable(CtrlVar,MUA,h); title(sprintf('h at t=%g',time))
    hold on    
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    %Plot_sbB(CtrlVar,MUA,s,b,B) ; title(sprintf('time=%g',time))
    
    
    subplot(4,1,2)
    QuiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),ub,vb,CtrlVar);
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    hold off
    
    subplot(4,1,3)
    PlotMeshScalarVariable(CtrlVar,MUA,dhdt);   title(sprintf('dhdt at t=%g',time))
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    
    subplot(4,1,4)
    PlotMeshScalarVariable(CtrlVar,MUA,ab);   title(sprintf('ab at t=%g',time))
    hold on
    
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    hold off
    
    
    x=MUA.coordinates(:,1);
    y=MUA.coordinates(:,2);
    
    Fb=scatteredInterpolant(x,y,b);
    Fs=Fb ; Fs.Values=s;
    
    xProfile=min(x):1000:max(x);
    
    yCentre=40e3+xProfile*0;
    sProfile=Fs(xProfile,yCentre);
    bProfile=Fb(xProfile,yCentre);
    
    BProfile=MismBed(xProfile,yCentre);
    
    
    fig200=figure(200);
    fig200.Position=[1200 50 figsWidth 2*figHeights];
    
    plot(xProfile/1000,sProfile,'b')
    hold on
    plot(xProfile/1000,bProfile,'b')
    plot(xProfile/1000,BProfile,'k')
    title(sprintf('t=%g',time))
    hold off
    
    
    fig300=figure(300);
    fig300.Position=[1200 700 figsWidth figHeights];
    PlotMuaMesh(CtrlVar,MUA)
    hold on 
    
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r','LineWidth',2);
    title(sprintf('t=%g',time))
    hold off
    
    drawnow
    %%
end


end