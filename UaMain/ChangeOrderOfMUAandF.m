



%%

function [UserVar,MUA,F,BCs]=ChangeOrderOfMUAandF(CtrlVar,UserVar,MUA,F,BCs,To)


%% Changes the order of the elements in the mesh.
%
%  To change all the elements of an existing mesh with 6-node elements to 3-node elements:
%
%    [UserVar,MUA,F,BCs]=ChangeOrderOfMUAandF(CtrlVar,UserVar,MUA,F,BCs,3);
%
%  To change all the elements of an existing mesh with 3-node elements to 10-node elements:
%
%    [UserVar,MUA,F,BCs]=ChangeOrderOfMUAandF(CtrlVar,UserVar,MUA,F,BCs,10);
%
% When reducing the order of elements, not interpolation is required as all nodes of the new mesh are contained within the
% original mesh.
%
% When increasing the order of elements, linear interpolation is used onto the new nodes using the MATLAB scatteredInterpolant.
%
%%

%%
% Klear ; load("TestSave.mat","CtrlVar","UserVar","MUA","F","BCs","l") ; [UserVar,MUA,F,BCs]=MeshTri6toTri3(CtrlVar,UserVar,MUA,F,BCs,3) ;



[Nele,nod]=size(MUA.connectivity);

if ~ismember(To,[3 6 10])
    fprintf('Changing to %-i nod element not possible\n',To)
    error('ChangeElementType')
end

From=nod;

if From==To ; return ; end

CtrlVar.CalcMUA_Derivatives=false;
CtrlVar.FindMUA_Boundary=true;
CtrlVar.MUA.MassMatrix=false ;
CtrlVar.MUA.StiffnessMatrix=false;
CtrlVar.MUA.DecomposeMassMatrix=false ;
CtrlVar.MUA.DecomposeMassMatrix=false ;
CtrlVar.Parallel.uvAssembly.spmd.isOn=false ;
CtrlVar.Parallel.uvhAssembly.spmd.isOn=false ;
CtrlVar.InfoLevelAdaptiveMeshing=0 ;

%%

FindOrCreateFigure("MUA in") ; PlotMuaMesh(CtrlVar,MUA)
UaPlots(CtrlVar,MUA,F,"-uv-",FigureTitle=" vel in ")

%%

fieldNames={"s";"b";"h";"S";"B";"ub";"vb";"ud";"vd";"AGlen";"n";"C";"m";"as";"ab";"dasdh";"dabdh";"rho";"dhdt";"dubdt";"dvbdt";"duddt";"dvddt"};


if (From==6 || From==10)  && To==3  % reducing the order

    % 6 to 3

    if From==6
        MUA.connectivity=MUA.connectivity(:,[1,3,5]);
    elseif From==10
        MUA.connectivity=MUA.connectivity(:,[1,4,7]);
    end
    Nodes=MUA.connectivity(:) ;
    C=unique(Nodes);
    B=zeros(MUA.Nnodes,1) ;
    B(C(1:numel(C)))=(1:numel(C));

    % for l=1:numel(C)
    %     B(C(l))=l;
    % end

    MUA.connectivity=B(MUA.connectivity) ;
    F.x=F.x(C);
    F.y=F.y(C);

    for k=1:numel(fieldNames)

        currentFieldName = fieldNames{k};
        F.(currentFieldName)=F.(currentFieldName)(C);

    end

    MUA.coordinates=[F.x(:) F.y(:)];

    CtrlVar.TriNodes=3;
    MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);

elseif  From==3 && (To==6 || To==10)    % increasing the order (I think this might also do the 6 to 10 node case)


    % One could do something clever using the structure of the mesh, but I suspect that I can't beat MATLAB scatteredinterpolant

    CtrlVar.TriNodes=To;
    [MUA.coordinates,MUA.connectivity]=ChangeElementType(MUA.coordinates,MUA.connectivity,CtrlVar.TriNodes) ;
    MUA=CreateMUA(CtrlVar,MUA.connectivity,MUA.coordinates);

    FInterpolant=scatteredInterpolant();
    FInterpolant.Points=[F.x F.y];

    F.x=MUA.coordinates(:,1) ; F.y=MUA.coordinates(:,2);
    for k=1:numel(fieldNames)

        currentFieldName = fieldNames{k};
        FInterpolant.Values=F.(currentFieldName);
        F.(currentFieldName)=FInterpolant(F.x,F.y);


    end


else

    error("Case not implemented")

end

[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
[UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F) ;

%%
FindOrCreateFigure("MUA out") ; PlotMuaMesh(CtrlVar,MUA)
UaPlots(CtrlVar,MUA,F,"-uv-",FigureTitle=" vel out ")


%%