%%
fid = fopen('GmshInitialGeoInputFile.msh');

for I=1:4
tline = fgetl(fid); 

end


Nodes = fscanf(fid, '%g',1);
NodeList= fscanf(fid, '%g',[4 Nodes]);

NodeList=NodeList';
tline = fgetl(fid);  disp(tline)
tline = fgetl(fid);  disp(tline)
tline = fgetl(fid);  disp(tline)


Ele = fscanf(fid, '%g',1);
EleList= fscanf(fid, '%g',[6 Ele]);

fclose(fid);
%%