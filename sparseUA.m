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



% if numel(varargin{1})>100000
%    
%     t2=tic ; S=sparse2(varargin{:}); t2=toc(t2) ; 
%     t=tic ; S=sparse(varargin{:}); t=toc(t) ; 
%     
%     fprintf(' Sparse2 in %f sec \t \t sparse in %f sec. \n',t2,t)
%     fprintf(' \n ')
% end

end