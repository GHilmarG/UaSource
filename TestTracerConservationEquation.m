
%%
UserVar.SUPG=1;  % True if using the SUPG
MainFigure=2;    % Just the number of the main figure windwo
DoPlots=1;

%% Main input variables
hStep=1; hmin=1;

kappa=0;     U0=10;  % No diffusion, infinite Peclet number
%kappa=1e-4;  U0=10;  % advection dominated, large Peclet number.
%kappa=1;  U0=10;  % a bit of both,
%kappa=1;     U0=0;   % no advection, Peclet number equal to zero

%% Create mesh
CtrlVar=Ua2D_DefaultParameters(); %
CtrlVar.TriNodes=3;
CtrlVar.MeshSizeMax=1;
CtrlVar.MeshSizeMin=0.1;
CtrlVar.MeshSize=0.1;
MeshBoundaryCoordinates=[-10 -1 ; -10 1 ; 10 1 ; 10 -1];
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;

[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotMuaMesh(CtrlVar,MUA)

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);   

%% Define velocities and initial tracer (thickness) distribution 
u0=zeros(MUA.Nnodes,1)+U0 ;  v0=zeros(MUA.Nnodes,1) ;
a0=zeros(MUA.Nnodes,1) ;
u1=u0 ; v1=v0 ; a1=a0;
c0=zeros(MUA.Nnodes,1) ;


c0(x<0)=hmin+hStep ;  c0(x>=0)=hmin ;  % Here I define a step change in thickness at x=0

figure(10) ; PlotMeshScalarVariable(CtrlVar,MUA,u0) ; 

% calculate CFL condition and make sure time step is small enough
Tarea=TriAreaFE(MUA.coordinates,MUA.connectivity);
Tlength=sqrt(2*Tarea) ;



switch CtrlVar.TriNodes

    case 6
        Tlength=Tlength/2; 
    case 10
        Tlength=Tlength/3;
end

CFL=min(Tlength)/max(u0);
dt=0.01; 
%dt=min([dt,CFL/2]);


time=0;

if ~UserVar.SUPG

    CtrlVar.SUPG.beta0=0 ; CtrlVar.SUPG.beta1=0 ;
else
    CtrlVar.SUPG.beta0=1/2 ; CtrlVar.SUPG.beta1=0 ;  % 
end

%% Define BCs
BCsTracer=BoundaryConditions;
I=intersect(find(MUA.coordinates == -10 ),MUA.Boundary.Nodes);
BCsTracer.hFixedNode=I; BCsTracer.hFixedValue=BCsTracer.hFixedNode*0+2; 




if DoPlots
    
    PlotBoundaryConditions(CtrlVar,MUA,BCsTracer); 
    
    Y0ind=abs(y)<min(Tlength)/10;
    Figh=figure(MainFigure) ;
    Figh.Position=[1200 800 1700 800];
    subplot(2,1,1)
    PlotMeshScalarVariable(CtrlVar,MUA,c0) ;
    subplot(2,1,2)
    plot(x(Y0ind),c0(Y0ind)) ;
end

speed=sqrt(u0.*u0+v0.*v0); 
fprintf(' Local element Peclet number |v| l/2 kappa = %f \n',max(speed)*min(Tlength)./2/mean(kappa)) 

%% And now finally calculate the trace evolution 
for I=1:25
    
    CtrlVar.theta=0.75;  % Note, it is quite possible that backward Euler is the best approach here
                      % So instead of using the default value of theta=1/2, set theta=1. 
                      % 
    [UserVar,c1,lambdah]=TracerConservationEquation(UserVar,CtrlVar,MUA,dt,c0,u0,v0,a0,u1,v1,a1,kappa,BCsTracer) ;
    
    
    time=time+dt;
    if DoPlots
        figure(MainFigure) ;
        subplot(2,1,1)
        PlotMeshScalarVariable(CtrlVar,MUA,c1) ;
        subplot(2,1,2)
        %plot(x(Y0ind),c1(Y0ind),'o-') ;
        plot(x,c1,'.') ;
        
        if UserVar.SUPG
            title(sprintf('c at t=%f with SUPG',time))
        else
            title(sprintf('c at t=%f without SUPG ',time))
        end
    end
    
    %     prompt = 'Continue? Y/N [Y]: ';
%     str = input(prompt,'s');
%     if isempty(str)
%         str = 'Y';
%     end
%     
%     if contains(lower(str),'n')
%         break
%     end
    
    c0=c1;
end

