function MESH = mtriangle(opts, tristr)
%MTRIANGLE calls the 'triangle' library by Jonathan Shewchuk to generate an
%unstructured triangular mesh. 
%
% Usage:
%  MESH = mtriangle(opts, triangle_input)
%
% Input:
%  opts is a structure containing options that influence triangle. Below
%  are the options and their default values
%
%   opts.element_type     = 'tri3' % use 'tri3', 'tri6', or 'tri7'
%   opts.triangulate_poly = 1;     % Triangulates a Planar Straight Line Graph
%   opts.gen_edges        = 0;     % 1: return edge information
%   opts.gen_neighbors    = 0;     % 1: return element neighbor information
%   opts.gen_elmarkers    = 0;     % 1: return element markers
%   opts.gen_boundaries   = 1;     % 1: return node markers
%   opts.min_angle        = 33;    % minimum triangle angle
%   opts.max_tri_area     = Inf;   % maximum triangle area
%   opts.ignore_holes     = 0;     % 
%   opts.exact_arithmetic = 1;     % 0: do not use Triangle's exact arithmetic
%   opts.zero_numbering   = 0;     % 1: use zero-based index numbering
%   opts.other_options    = '';    % other triangle options
% 
%  triangle_input is a structure defining the input data passed to
%  triangle. For a complete explanation read the triangle documentation.
%    
%   tristr.points
%   tristr.segments
%   tristr.segmentmarkers
%   tristr.regions
%   tristr.pointmarkers
%   tristr.pointattributes
%   tristr.holes
%   tristr.triangles
%   tristr.triangleattributes
%   tristr.trianglearea
%
% Example:
%
%  points   = [0 0; 1 0; 1 1; 0 1]';
%  segments = [1 2; 2 3; 3 4; 4 1]';
%  opts.element_type     = 'tri3';
%  opts.min_angle        = 30;
%  opts.max_tri_area     = 0.001;
%  tristr.points         = points;
%  tristr.segments       = int32(segments);
%  MESH = mtriangle(opts, tristr);

% Copyright 2012, Marcin Krotkiewski, University of Oslo

%% parameter checks

if ~isfield(tristr, 'points')
    tristr.points = [];
end
if ~isfield(tristr, 'segments')
    tristr.segments = [];
end
if ~isfield(tristr, 'segmentmarkers')
    tristr.segmentmarkers = [];
end
if ~isfield(tristr, 'regions')
    tristr.regions = [];
end
if ~isfield(tristr, 'pointmarkers')
    tristr.pointmarkers = [];
end
if ~isfield(tristr, 'pointattributes')
    tristr.pointattributes = [];
end
if ~isfield(tristr, 'holes')
    tristr.holes=[];
end
if ~isfield(tristr, 'triangles')
    tristr.triangles=int32([]);
end
if ~isfield(tristr, 'triangleattributes')
    tristr.triangleattributes=[];
end
if ~isfield(tristr, 'trianglearea')
    tristr.trianglearea=[];
end

opts = analyze_options(opts);
tri_flag = generate_options_string(opts);

%% run triangle
if exist(['triangle.' mexext]) == 3
    %% use the mex file
    [MESH.NODES, MESH.ELEMS, MESH.elem_markers, MESH.node_markers, ...
        MESH.EDGES, MESH.edge_markers, ...
        MESH.SEGMENTS, MESH.segment_markers, ...
        MESH.NEIGHBORS] = ...
        triangle(tri_flag, tristr.points, tristr.pointmarkers, tristr.pointattributes, ...
        tristr.segments, tristr.segmentmarkers, tristr.holes, tristr.regions, ...
        tristr.triangles, tristr.triangleattributes, tristr.trianglearea);   
else
    % only mex file in standalone distribution
    if triangle_standalone
        error('MEX function not foun');
    end
    
    %% triangle executable
    %  first find triangle in the path
    if exist('triangle') == 2
        texec = 'triangle';
        print_message('Using triangle found in the path.\n');
    else
        
        milamin_data = getappdata(0, 'milamin_data');
        if ispc
            texec = [milamin_data.path 'ext/triangle/triangle.exe'];
        else
            texec = [milamin_data.path 'ext/triangle/triangle.' lower(computer)];
        end
        if ~exist(texec)
            error([mfilename ': triangle executable not found: ' texec]);
        end
    end
    
    % delete old files
    delete([pwd filesep opts.model_name '.*.poly']); 
    delete([pwd filesep opts.model_name '.*.edge']); 
    delete([pwd filesep opts.model_name '.*.ele']); 
    delete([pwd filesep opts.model_name '.*.node']); 
    delete([pwd filesep opts.model_name '.*.neigh']); 
    
    % write input files
    triangle_write([pwd filesep opts.model_name], tristr.points, [tristr.segments; tristr.segmentmarkers], tristr.regions);
    
    [status,result] = system([texec ' -' tri_flag ' ' opts.model_name '.poly']);
    if status
        error([mfilename ': triangle executable failed: ' result]);
    end
    [MESH.NODES, MESH.ELEMS, MESH.elem_markers, MESH.node_markers, ...
        MESH.EDGES, MESH.edge_markers, ...
        MESH.SEGMENTS, MESH.segment_markers, ...
        MESH.NEIGHBORS] = ...
        triangle_read([pwd filesep opts.model_name '.1']);
end

if isempty(MESH.elem_markers)
    MESH.elem_markers = zeros(1, length(MESH.ELEMS));
end

if isfield(MESH, 'EDGES') & isempty(MESH.EDGES)
    MESH = rmfield(MESH, 'EDGES');
end
if isfield(MESH, 'edge_markers') & isempty(MESH.edge_markers)
    MESH = rmfield(MESH, 'edge_markers');
end
if isfield(MESH, 'NEIGHBORS') 
    if isempty(MESH.NEIGHBORS)
        MESH = rmfield(MESH, 'NEIGHBORS');
    else
        temp = MESH.NEIGHBORS==-1;
        if ~opts.zero_numbering
            MESH.NEIGHBORS(temp)=0;
        end
    end
end

% extra functionality comes with MILAMIN_v2
if ~triangle_standalone
    MESH = mesh_info(MESH,0);

    % - triangle does generate edge information, but does not
    %   create ELEMS_EDGES, i.e., edge-based element definitions
    % - triangle returns only two-node edges and disregards midedge nodes
    %   for higher-order elements.
    el_info = element_info(opts.element_type);
    if el_info.order>1 | opts.gen_edges
        MESH = mesh_find_edges(MESH);
    end
    
    if strcmp(opts.element_type, 'tri7')
        MESH = mesh_convert(MESH, 7);
    end
else
  if strcmp(opts.element_type, 'tri7')
    nel = size(MESH.ELEMS, 2);
    ncorners = 3;
    MESH.ELEMS(end+1,:)  = max(MESH.ELEMS(:))+int32([1:nel]);
    MESH.NODES = [MESH.NODES [...
			      mean(reshape(MESH.NODES(1, MESH.ELEMS(1:ncorners,:)), ncorners, nel));...
			      mean(reshape(MESH.NODES(2, MESH.ELEMS(1:ncorners,:)), ncorners, nel))]];
    if isfield(MESH, 'node_markers')
      MESH.node_markers = [MESH.node_markers zeros(1,nel)];
    end
  end
end

if ~opts.gen_edges & isfield(MESH, 'EDGES')
    MESH = rmfield(MESH, {'EDGES' 'edge_markers' 'ELEMS_EDGES'});
   if ~triangle_standalone
        MESH = mesh_info(MESH,0);
   end
end

%% parameter analysis functions
    function opts = analyze_options(opts)
        if ~isfield(opts, 'model_name')
            opts.model_name = 'model';
        end
        if ~isfield(opts, 'triangulate_poly')
            opts.triangulate_poly = 1;
        end
        if ~isfield(opts, 'element_type')
            opts.element_type = 'tri3';
        end
        if ~isfield(opts, 'gen_edges')
            opts.gen_edges        = 0;
        end
        if ~isfield(opts, 'gen_neighbors')
            opts.gen_neighbors    = 0;
        end
        if ~isfield(opts, 'gen_elmarkers')
            opts.gen_elmarkers    = 0;
        end
        if ~isfield(opts, 'gen_boundaries')
            opts.gen_boundaries   = 1;
        end
        if ~isfield(opts, 'min_angle')
            opts.min_angle        = 33;
        end
        if ~isfield(opts, 'max_tri_area')
            opts.max_tri_area     = 0;
        end
        if ~isfield(opts, 'ignore_holes')
            opts.ignore_holes     = 0;
        end
        if ~isfield(opts, 'exact_arithmetic')
            opts.exact_arithmetic = 1;
        end
        if ~isfield(opts, 'zero_numbering')
            opts.zero_numbering   = 0;
        end
        if ~isfield(opts, 'other_options')
            opts.other_options    = '';
        end
    end 
    function opts_str = generate_options_string(opts)
        
        opts_str = [opts.other_options];
        if opts.triangulate_poly
            opts_str = [opts_str 'p'];
        end
        switch opts.element_type
            case 'tri3'
                opts_str = [opts_str 'o1'];
            case {'tri6' 'tri7'}
                opts_str = [opts_str 'o2'];
            otherwise
                error([mfilename ': unknown triangle element type.']);
        end
        if opts.gen_edges
            opts_str = [opts_str 'e'];
        end
        if opts.gen_neighbors
            opts_str = [opts_str 'n'];
        end
        if opts.gen_elmarkers
            opts_str = [opts_str 'A'];
        end
        if ~opts.gen_boundaries
            opts_str = [opts_str 'B'];
        end
        opts_str = [opts_str 'q' num2str(opts.min_angle, '%.16f')];
        if opts.max_tri_area
            opts_str = [opts_str 'a' num2str(opts.max_tri_area, '%.16f')];
        end
        if opts.ignore_holes
            opts_str = [opts_str 'O'];
        end
        if ~opts.exact_arithmetic
            opts_str = [opts_str 'X'];
        end
        if opts.zero_numbering
            opts_str = [opts_str 'z'];
        end
    end
    function retval = triangle_standalone
        if isempty(getappdata(0, 'milamin_data'))
            retval = 1;
        else
            retval = 0;
        end
    end
end
