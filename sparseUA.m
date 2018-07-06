function S=sparseUA(varargin)

persistent  isSparse2



if isempty(isSparse2)
    if exist('sparse2','file')==3
        isSparse2=1 ;
    else
        isSparse2=0 ;
        
        fprintf('INFO: The function sparse2 is not known to your matlab installation.\n')
        fprintf('Matlab in-built sparse function used instead.\n')
        fprintf(' Consider installing SuiteSparse as this might speed up the runs. \n')
        fprintf(' see: http://faculty.cse.tamu.edu/davis/suitesparse.html\n')
    end
end


if isSparse2
    S=sparse2(varargin{:});
else
    S=sparse(varargin{:});
end


end