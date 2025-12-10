function Lima=ReadLima(FileName)

% Lima=ReadLima(FileName)
% To plot:
% figure ; image(Lima.x,Lima.y,Lima.Image) ; axis xy
%%
CurDir=pwd;

if nargin<1
    AntarcticGlobalDataSets=getenv('AntarcticGlobalDataSets');
    cd(AntarcticGlobalDataSets)
    fprintf('Reading Lima Composite')
    load AntarcticLimaComposite.mat
    fprintf(' done \n ')
else
    fprintf('Reading Lima Composite %s',FileName)
    load(FileName)
    fprintf(' done \n ')
end

cd(CurDir)

Lima.Image=AntarticaLimaComposite;
Lima.R=R;
Lima.bbox=bbox;
Lima.x=linspace(bbox(1,1),bbox(2,1),size(AntarticaLimaComposite,1));
Lima.y=linspace(bbox(2,2),bbox(1,2),size(AntarticaLimaComposite,2));

end

