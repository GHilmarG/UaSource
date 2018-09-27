function varargout=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,varargin)

%%
% varargout=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,varargin)
%
% OutsideValues not used.  (This was needed in older versions of Matlab before
% scatteredinterpolant started to support nearest-neighbour extrapolation)
%
% example:
% [h2,s2,rho2]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,h1,s1,rho1)
% interpolates h,s and rho, from FE Mesh MUA1 onto (x2,y2).

%%

nVarargsIn = length(varargin);

%fprintf('MapNodalVariablesFromMesh1ToMesh2 : Total number of varargin inputs = %d\n',nVarargsIn);

nVar=nVarargsIn;

varargout=cell(nVar,1);

if isempty(x2)
    for iOut=1:nVar
        varargout{iOut}=NaN;
    end
    return
end

% 
% if isempty(OutsideValues)
%     OutsideValues=zeros(nVar,1)+NaN;
% elseif numel(OutsideValues)~=nVar
%     if numel(OutsideValues)==1
%         OutsideValues=zeros(nVar,1);
%     else
%         fprintf('in MapNodalVariablesFromMesh1ToMesh2 the number elements in OutsideValues (%-i) must be equal the number of varargin cells (%-i) \n',numel(OutsideValues),nVar)
%         error('MapNodalVariablesFromMesh1ToMesh2:WrongNumberOfOutsideValues','Wrong number of elements in the vector OutSideValues')
%     end
% end

F = scatteredInterpolant();
F.Points=MUA1.coordinates;

F.Method='natural';
F.ExtrapolationMethod='nearest';

% The underlying triangulation is only done once, and when interpolating different fields
% only the values are changed.

%tri=TriFE(MUA1.connectivity);  % Only corner points
tri=CreateFEmeshTriRep(MUA1.connectivity,MUA1.coordinates);


tol = eps*1000;
[ID,d] = nearestNeighbor(tri,[x2(:) y2(:)]);
same=d<tol ;


for iVar=1:nVar
    
    if isempty(varargin{iVar})
        varargout{iVar}=[];
    else
        F.Values=double(varargin{iVar});
        
        varargout{iVar}=srMAP(x2,y2,F,ID,same);
        
        %varargout{iVar}=F(x2,y2);
    end
end

%% This is no longer needed because scatteredInterpolant now supports extrapolation
% %% If F.ExtrapolationMethod='none' then values outside of the convex hull are automatically
% % set to NaN. In many cases this is fine.  If, for example, the boundary and the topology of the FE mesh does
% % not change with time, there are no points outside of the mesh and this extrapolation option is irrelevant
% %
% % But in the general complex case where the FE mesh evolves with time one has to
% % locate all points outside of the boundary as defined in MUA1 (i.e. old mesh)
% % and set values outside to `OutsideValue'
% 
% if strcmp(F.ExtrapolationMethod,'none')
%     
%     if CtrlVar.Mesh1To2CheckForPointsOutsideMesh1AndInsideConvexHull
%         
%         % now check if there are locations inside of the convex hull but outside of the
%         % FE mesh boundaries.
%         % p=([x2(~I)  y2(~I)]) ; % only need to check for points inside of the convex hull
%         
%         p=([x2(:) y2(:)]);
%         node=MUA1.coordinates;
%         edge=MUA1.Boundary.Edges(:,[1 end]);  % here I assume straight element edges!
%         
%         [in,on]=inpoly(p,node,edge,CtrlVar.InpolyTol);
%         
%         % set all values outside of the FE boundary to the OutsideValue
%         for iVar=1:nVar
%             varargout{iVar}(~in)=OutsideValues(iVar);   % outside of mesh boundary
%         end
%         
%     end
%     
% end

%% set values where (x,y) is NaN to NaN
Inan=isnan(x2) | isnan(y2) ;
for iVar=1:nVar
    varargout{iVar}(Inan)=NaN;
end







end