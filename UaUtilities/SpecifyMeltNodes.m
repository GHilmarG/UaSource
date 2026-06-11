function [MeltNodes,GLgeo,GLnodes,GLele,OceanNodes,LakeNodes]=SpecifyMeltNodes(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele)

%%
%
%   [MeltNodes,NotMeltNodes,GLgeo,GLnodes,GLele,OceanNodes,LakeNodes]=SpecifyMeltNodes(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele)
%
% tries to come up with a resonable definition of nodes to which ocean-induced
% melt should be applied to

if nargin<4 || isempty(GLgeo) || isempty(GLnodes) || isempty(GLele)
    [GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
end

[OceanNodes,LakeNodes]=LakeOrOcean(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele);



if MUA.nod ==3
    
    % All ocean nodes, and only ocean nodes, are melt nodes
    MeltNodes=OceanNodes;
    
else
    
    switch lower(CtrlVar.MeltNodesDefinition)
        
        case 'node-wise'
            
            MeltNodes=OceanNodes;
            
            
            
        case 'edge-wise'
            
            % All ocean nodes are melt nodes.
            % In additon: If both corner nodes of an edge are melt nodes, then the nodes along the edge are as well.
            
            switch MUA.nod
                case 6
                    
                    I=AllElementsContainingGivenNodes(MUA.connectivity,OceanNodes);
                    
                    e=zeros(3,2);
                    e(1,:)=[1 3]; e(2,:)=[3 5] ; e(3,:)=[5 1];
                    
                    for ie=1:3
                        Iedge{ie}=all(GF.node(MUA.connectivity(:,e(ie,:)))<CtrlVar.GLthreshold,2);
                        Iedge{ie}=Iedge{ie} & I ;
                        Iedge{ie}=MUA.connectivity(Iedge{ie},2*ie);
                    end
                    
                    MeltNodes=unique([Iedge{1}(:);Iedge{2}(:);Iedge{3}(:);OceanNodes(:)]);
                    
                case 10
                    error('10 node case not yet implemented')
            end
            
        case 'element-wise'
            
            % All ocean nodes are melt nodes.
            % In additon: If all corner nodes are melt nodes, then all others are as well.
            
            I=AllElementsContainingGivenNodes(MUA.connectivity,OceanNodes);
            
            switch MUA.nod
                case 6
                    Iele=all(GF.node(MUA.connectivity(:,[1 3 5]))<CtrlVar.GLthreshold,2);
                case 10
                    Iele=all(GF.node(MUA.connectivity(:,[1 4 7]))<CtrlVar.GLthreshold,2);
            end
            Iele=Iele  & I ;  % ele containing ocean nodes with all corner nodes afloat
            
            NodesOfMeltEle=MUA.connectivity(Iele,:);  %
            MeltNodes=unique([NodesOfMeltEle(:);OceanNodes(:)]);
            
        otherwise
            error('DefineMeltNodes: Case not reckognised\m')
    end
end



%%
if CtrlVar.PlotMeltNodes && CtrlVar.doplots
    %%
    [GLgeo,GLinfo,GLnodes]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
    figure ; PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
    hold on
    x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);
    
    
    
    plot(x(OceanNodes)/CtrlVar.PlotXYscale,y(OceanNodes)/CtrlVar.PlotXYscale,'xb')
    I=GF.node<CtrlVar.GLthreshold; plot(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,'+r')
    plot(x(MeltNodes)/CtrlVar.PlotXYscale,y(MeltNodes)/CtrlVar.PlotXYscale,'ok')
    
    plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'g','LineWidth',1);
    xlabel(CtrlVar.PlotsXaxisLabel) ; ylabel(CtrlVar.PlotsYaxisLabel) ;
    title(['Melt nodes: time=',num2str(CtrlVar.time)])
    legend('Mesh','Ocean','afloat','melt nodes')
    %%
end
%%
end

