function Data=CompareRuns(files,Variable)

%%
%
% Plots differences between results in data files
%
% files     :  a cell array with files names
%
%
% Example  
%
%   files{1}="ResultsfileOne.mat"; 
%   files{2}="ResultsfileTwo.mat"; 
%   CompareRuns(files,'ubvb') ;
%
%
% Note: This is a very simple m-file with limited options. 
%
%%



N=numel(files);  Data=cell(N,1) ;

if ischar(files{1}) || isstring(files{1})
    
    
    for I=1:N
        Data{I}=load(files{I});
    end
    
    for I=1:N
        if ~isfield(Data{I},'CtrlVar')  && isfield(Data{I},'CtrlVarInRestartFile')
            Data{I}.CtrlVar=Data{I}.CtrlVarInRestartFile;
        end
    end
    
   
  
else
    fprintf("input needs to contain file names. \n")
    return
end

for I=1:N
    [cbar,xGL,yGL,xCF,yCF]=UaPlots(Data{I}.CtrlVar,Data{I}.MUA,Data{I}.F,Variable,FigureTitle=files{I}) ; 
end

for I=2:N
    if Data{1}.MUA.Nnodes==Data{I}.MUA.Nnodes

        RefData=Data{1} ; % This data in the first data file is used as a reference

        RefData.F.ub=Data{I}.F.ub-Data{1}.F.ub;
        RefData.F.vb=Data{I}.F.vb-Data{1}.F.vb;
        RefData.F.s=Data{I}.F.s-Data{1}.F.s;
        RefData.F.b=Data{I}.F.b-Data{1}.F.b;
        RefData.F.h=Data{I}.F.h-Data{1}.F.h;
        RefData.F.dhdt=Data{I}.F.dhdt-Data{1}.F.dhdt;

        [cbar,xGL,yGL,xCF,yCF]=UaPlots(RefData.CtrlVar,RefData.MUA,RefData.F,Variable,FigureTitle=files{I}+" minus "+files{1}) ;


    end
end

end
