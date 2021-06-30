function hTri=PlotFEmesh(coordinates,connectivity,CtrlVar,ElementList,varargin)
 
%
% PlotFEmesh(coordinates,connectivity)
% PlotFEmesh(coordinates,connectivity,CtrlVar)
% PlotFEmesh(coordinates,connectivity,CtrlVar,ElementList)   only plots elements in the list ind (ind is not a locical index)
%
% plots FE mesh
%
%       Default values:
%         CtrlVar.PlotLabels=0;
%         CtrlVar.MeshColor='k';
%         CtrlVar.NodeColor='k';
%         CtrlVar.PlotXYscale=1;
%         CtrlVar.PlotNodesSymbolSize=3;
%         CtrlVar.PlotNodesSymbol='o';
%         CtrlVar.PlotNodes=0;
%         CtrlVar.PlotMesh=1;
%         CtrlVar.time=NaN;
%         CtrlVar.PlotFEmeshAndSaveMesh=0;



persistent iCounter



if nargin>3 && ~isempty(ElementList)
    connectivity=connectivity(ElementList,:);
    if islogical(ElementList)
        ElementNumbers=find(ElementList);
    else
        ElementNumbers=ElementList;
    end
else
    ElementNumbers=1:size(connectivity,1);
end

ElementNumbers=ElementNumbers(:);
[Nele,nod]=size(connectivity);

if Nele==0
    hTri=[];
    return
end

if nargin < 3
    CtrlVar.PlotLabels=0;
    CtrlVar.MeshColor='k';
    CtrlVar.NodeColor='k';
    CtrlVar.PlotXYscale=1;
    CtrlVar.PlotNodesSymbolSize=3;
    CtrlVar.PlotNodesSymbol='o';
    CtrlVar.PlotNodes=1;
    CtrlVar.time=NaN;
    CtrlVar.FEmeshPlotTitle=[];
    CtrlVar.PlotFEmeshAndSaveMesh=0;
    CtrlVar.PlotsXaxisLabel='x';
    CtrlVar.PlotsYaxisLabel='y';
else
    
    if ~isfield(CtrlVar,'PlotLabels') ; CtrlVar.PlotLabels=0; end
    if ~isfield(CtrlVar,'MeshColor') ; CtrlVar.MeshColor='k'; end
    if ~isfield(CtrlVar,'NodeColor') ; CtrlVar.NodeColor='k'; end
    if ~isfield(CtrlVar,'PlotXYscale') ; CtrlVar.PlotXYscale=1; end
    if ~isfield(CtrlVar,'PlotNodesSymbol') ; CtrlVar.PlotNodesSymbol='o'; end
    if ~isfield(CtrlVar,'PlotNodesSymbolSize') ; CtrlVar.PlotNodesSymbolSize=3; end
    if ~isfield(CtrlVar,'PlotNodes') ; CtrlVar.PlotNodes=1; end
    if ~isfield(CtrlVar,'time') ; CtrlVar.time=NaN; end
    if ~isfield(CtrlVar,'FEmeshPlotTitle') ; CtrlVar.FEmeshPlotTitle=[];end
    if ~isfield(CtrlVar,'PlotFEmeshAndSaveMesh') ; CtrlVar.PlotFEmeshAndSaveMesh=0;end
    if ~isfield(CtrlVar,'PlotsXaxisLabel') ; CtrlVar.PlotsXaxisLabel='x'; end
    if ~isfield(CtrlVar,'PlotsYaxisLabel') ; CtrlVar.PlotsYaxisLabel='x'; end
    
    
end


if ~isfield(CtrlVar,'PlotEleLabels')
    if CtrlVar.PlotLabels
        CtrlVar.PlotEleLabels=1;
    else
        CtrlVar.PlotEleLabels=0;
    end
end


if ~isfield(CtrlVar,'PlotNodalLabels')
    if CtrlVar.PlotLabels
        CtrlVar.PlotNodalLabels=1;
    else
        CtrlVar.PlotNodalLabels=0;
    end
end


if isempty(iCounter) ; iCounter=0 ; end

if CtrlVar.PlotFEmeshAndSaveMesh
     if iCounter==0 && exist('ResultsFiles','dir')~=7 
        mkdir('ResultsFiles') ;
    end
    iCounter=iCounter+1;
    FileName=['ResultsFiles/',sprintf('%07i',iCounter),'-Mesh-',CtrlVar.Experiment];
    fprintf('Saving mesh in %s \n',FileName)
    save(FileName,'coordinates','connectivity')
end


coordinates=coordinates/CtrlVar.PlotXYscale;

%FEmeshCPT=CreateFEmeshCornerPointTriangulation(connectivity,coordinates);
TR=CreateFEmeshTriRep(connectivity,coordinates);

if numel(varargin)==1
   varargin={'color',varargin{1}};
end

hTri=triplot(TR,'color',CtrlVar.MeshColor,varargin{:}) ;
hold on

if nargin>3
    nodes=unique(connectivity);
else
    nodes=1:length(coordinates);
end

Nnodes=numel(nodes);


if CtrlVar.PlotNodes==1
    plot(coordinates(nodes,1),coordinates(nodes,2),CtrlVar.PlotNodesSymbol,...
        'color',CtrlVar.NodeColor,'MarkerSize',CtrlVar.PlotNodesSymbolSize)
end


if CtrlVar.PlotNodalLabels
    fprintf(' labeling nodes\n')
    vxlabels = arrayfun(@(n) {sprintf(' %d', n)}, nodes(:));
    text(coordinates(nodes,1),coordinates(nodes,2),vxlabels, 'FontWeight', 'bold', 'HorizontalAlignment','center', 'BackgroundColor', 'none','color','blue');
end
%numtri = size(FEmeshTriRep,1);


if CtrlVar.PlotEleLabels
    fprintf(' labeling elements\n')
    ic = incenter(TR);
    trilabels = arrayfun(@(x) {sprintf('%d', x)}, ElementNumbers);
    text(ic(:,1), ic(:,2), trilabels, 'FontWeight', 'bold','HorizontalAlignment', 'center','Color', 'red');
end


[Nele,nod]=size(connectivity);

if isempty(CtrlVar.FEmeshPlotTitle)
    if isnan(CtrlVar.time)
        title(sprintf(' Ele=%-i,  Nodes=%-i,  nod=%-i',Nele,Nnodes,nod))
    else
        title(sprintf('time=%-g \t  Ele=%-i,  Nodes=%-i,  nod=%-i',CtrlVar.time,Nele,Nnodes,nod))
    end
else    
    title(sprintf('time=%-g \t  Ele=%-i,  Nodes=%-i,  nod=%-i',CtrlVar.time,Nele,Nnodes,nod))
end

xlabel(CtrlVar.PlotsXaxisLabel) ; ylabel(CtrlVar.PlotsYaxisLabel)
axis equal

end




