



filename='Ex3a3D-StraightChannelWidth40GLposition.mat';
filename='Ex3a3D-StraightChannelWidth50GLposition.mat';
%filename='Ex3a3D-StraightChannelWidth60GLposition.mat';
%filename='Ex3a3D-StraightChannelWidth70GLposition.mat';
%filename='GLpositionEx3aParameterSetOne.mat';
filename='GLpositionEx3aParameterSetTwo.mat';
filename='M3DinitialSteadyStateGLposition.mat';
%filename='GLposition.mat';
load(filename,'GLposMax','GLposMin','GLposMean','GLstd','GLtime') 

H0=1000; T0=1;

figure 
plot(GLtime/T0,GLposMax/1000,'r.')
hold on
plot(GLtime/T0,GLposMin/1000,'b.')
plot(GLtime/T0,GLposMean/1000,'g')
legend('GLmax','GLmin','GLmean','Location','NorthWest')
ax=axis;


