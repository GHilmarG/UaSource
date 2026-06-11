function [Lon,Lat,BasID] = ReadDrainageBasinsZwally

%% Initialize variables.
%filename = '/data/dataphy/janryd69/Antarctic_datasets/Antarctica_DrainageBasins_Zwally/Ant_Full_DrainageSystem_Polygons.txt';
%C:\cygwin64\home\Hilmar\ghg\Ua\Antarctic Global Data Sets\Antarctic Basins\Antarctica_DrainageBasins_Zwally

filename = 'Ant_Full_DrainageSystem_Polygons.txt';
startRow = 8;

%% Format string for each line of text:
formatSpec = '%15f%15f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
Lat = dataArray{:, 1};
Lon = dataArray{:, 2};
BasID = dataArray{:, 3};

%% Clear temporary variables
clearvars filename startRow formatSpec fileID dataArray ans;