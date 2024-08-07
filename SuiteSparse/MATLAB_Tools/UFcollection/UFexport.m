function UFexport (list, check, tmp)
%UFEXPORT export to Matrix Market and Rutherford/Boeing formats
%
% Example:
%   UFexport ;          % export the entire collection
%   UFexport (list) ;   % just export matrices whose id's are given in the list
%   UFexport (list, 'check') ;      % also read them back in, to check
%
% If the list is empty, all matrices in the collection are exported.
% A 3rd argument tmp changes the tmp directory for UFread.
%
% See also UFget, UFwrite, RBio, mwrite.


% Copyright 2006-2007, Timothy A. Davis

%-------------------------------------------------------------------------------
% get the input arguments
%-------------------------------------------------------------------------------

index = UFget ;
nmat = length (index.nrows) ;

if (nargin < 1 || isempty (list))
    list = 1:nmat ;
end

check = ((nargin > 1) && strcmp (check, 'check')) ;

if (nargin < 3)
    tmp = '' ;
end

%-------------------------------------------------------------------------------
% determine the top-level directory to use
%-------------------------------------------------------------------------------

[url topdir] = UFlocation ;
fprintf ('\nExport to topdir: %s\ncheck: %d\ntmp: %s\n', topdir, check, tmp) ;
for id = list
    fprintf ('%4d : %s/%s\n', id, index.Group {id}, index.Name {id}) ;
end

%-------------------------------------------------------------------------------
% export the matrices
%-------------------------------------------------------------------------------

format = { 'MM' , 'RB' } ;

for id = list

    % get the MATLAB version
    clear Problem
    Problem = UFget (id, index) ;
    disp (Problem) ;

    % create the MM and RB versions
    for k = 1:2
        fprintf ('Exporting to %s format ...\n', format {k}) ;
        if (nnz (Problem.A) < 1e8)
            UFwrite (Problem, [topdir format{k}], format{k}, 'tar') ;
        else
            % the MATLAB tar has problems with huge files
            fprintf ('File to big for MATLAB tar\n') ;
            UFwrite (Problem, [topdir format{k}], format{k}) ;
        end
    end

    % check the new MM and RB versions
    if (check)
        for k = 1:2
            fprintf ('Reading %s format ...\n', format {k}) ;
	    try
		if (isempty (tmp))
		    P2 = UFread ([topdir format{k} '/' Problem.name]) ;
		else
		    P2 = UFread ([topdir format{k} '/' Problem.name], tmp) ;
		end
	    catch
		% The Problem may be too large for two copies to be in the
		% MATLAB workspace at the same time.  This is not an error,
		% but it means that the Problem cannot be checked.
		P2 = [ ] ;
		fprintf ('Unable to read %s/%s\n', format {k}, Problem.name) ;
		fprintf ('%s\n', lasterr) ;
	    end
            fprintf ('Comparing MATLAB and %s format ...\n', format {k}) ;
	    if (~isempty (P2) && ~isequal (Problem, P2))
                Problem
                P2
		warning ('%s version mismatch: %s\n', format {k}, Problem.name);
                e = norm (Problem.A - P2.A, 1) ;
                a = norm (Problem.A,1) ;
                fprintf ('norm (A1-A2,1): %g  relative: %g\n', e, e/a) ;
	    end
	    clear P2
	end
        fprintf ('OK.\n') ;
    end
end

