function varargout=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,MUA1,x2,y2,varargin)

% varargout=srMapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,varargin)
%
% example:
% [h2,s2,rho2]=srMapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,h1,s1,rho1)
% interpolates h,s and rho, from FE Mesh MUA1 onto (x2,y2).

nVar = length(varargin);

varargout=cell(nVar,1);

tri=CreateFEmeshTriRep(MUA1.connectivity,MUA1.coordinates);

%% step 1: find which nodes are new in mesh2

tol = 1e-6;
[IDsame,d] = nearestNeighbor(tri,[x2(:) y2(:)]);
same=d<tol ;

nmap = sum(~same); % number of nodes that need to be mapped

if nmap==0 % if the nodes are all in the same locations, map directly and exit
    % this should work even if the node numbering has changed
    for iVar=1:nVar
        
        if isempty(varargin{iVar})
            varargout{iVar}=[];
        else
            varargout{iVar} = zeros(length(x2),1);
            varargout{iVar}(same) = varargin{iVar}(IDsame(same));
        end
    end
    return;
end


%% step 2: for new nodes - find their parent element and barycentric coords

xloc = x2(~same);
yloc = y2(~same);

% NB: the triangulation is reduced to 3 node no matter what the input is
% but this shouldn't matter


Outside = false(nmap,1);

% each row of ID is the element ID of the triangle enclosing that point
% each row of B are the barycentric coordinates of that point within
% element number ID
[ID,B] = pointLocation(tri,[xloc yloc]);

% PROBLEM!! many new nodes are along edges... these generate a NaN value in
% the call to pointLocation because they are in two elements at once.
% need to somehow move the point arbirtary small distance to get it
% definitively inside a triangle
% this should only affect the solution by a very tiny amount as long as the
% element size >> tol (micro meters)

if any(isnan(ID))
    % first try moving up right
    [tempID,tempB] = pointLocation(tri,[xloc(isnan(ID))+tol yloc(isnan(ID))+tol]);
    B(isnan(ID),:) = tempB;
    ID(isnan(ID)) = tempID;
    if any(isnan(ID))
        % next try moving down left
        [tempID,tempB] = pointLocation(tri,[xloc(isnan(ID))-tol yloc(isnan(ID))-tol]);
        B(isnan(ID),:) = tempB;
        ID(isnan(ID)) = tempID;
        if any(isnan(ID))
            % next try moving up left
            [tempID,tempB] = pointLocation(tri,[xloc(isnan(ID))-tol yloc(isnan(ID))+tol]);
            B(isnan(ID),:) = tempB;
            ID(isnan(ID)) = tempID;
            if any(isnan(ID))
                % finally try moving down right
                [tempID,tempB] = pointLocation(tri,[xloc(isnan(ID))+tol yloc(isnan(ID))-tol]);
                B(isnan(ID),:) = tempB;
                ID(isnan(ID)) = tempID;
                if any(isnan(ID))
                    % this point must lie outside of the triangulation
                    % altogether so mark it to deal with for later
                    Outside(isnan(ID)) = true;
                    ID(isnan(ID)) = 1; % arbitrary - will be fixed later
                end
            end
        end
    end
end


%% step 3: get the shape functions at each new node


newvals = zeros(nmap,1);

sfun = sr_shape_fun(B,CtrlVar.TriNodes);

%% step 4: use the shape functions to map new nodal values for each argin
for iVar=1:nVar
    
    if isempty(varargin{iVar})
        varargout{iVar}=[];
    else
        
        Fnode = reshape(varargin{iVar}(MUA1.connectivity,1),MUA1.Nele,MUA1.nod);
        
        for ii = 1:nmap
            newvals(ii) = Fnode(ID(ii),:)*sfun(ii,:)';
        end
        
        temp = IDsame(~same);
        newvals(Outside) = mean(varargin{iVar}(temp(Outside),:),2);
        varargout{iVar} = zeros(length(x2),1);
        varargout{iVar}(same) = varargin{iVar}(IDsame(same));
        varargout{iVar}(~same) = newvals;
        
    end
end
