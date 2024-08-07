function C = cs_transpose (A)                                               %#ok
%CS_TRANSPOSE transpose a real sparse matrix.
%   C = cs_transpose(A), computes C = A' where A must be sparse and real.
%
%   Example:
%       Prob = UFget ('HB/ibm32') ; A = Prob.A ;
%       C = cs_transpose (A) ;
%       C-A'
%
%   See also TRANSPOSE, CTRANSPOSE.

%   Copyright 2006-2007, Timothy A. Davis.
%   http://www.cise.ufl.edu/research/sparse

error ('cs_transpose mexFunction not found') ;


