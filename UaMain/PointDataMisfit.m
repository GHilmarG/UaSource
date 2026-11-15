function [J,dJdf,d2Jdf2]=PointDataMisfit(O,d,dErr,f)

%%
%
%   [J,dJdf,d2Jdf2]=PointDataMisfit(O,d,dErr,f)
%
% Misfit between point measurements, d, and a nodal finite-element field, f, evaluated at the measurement locations, with
% uncorrelated measurement errors:
%
% $$ J = \frac{1}{2} \, (\mathbf{O}\mathbf{f}-\mathbf{d})^T \, \mathbf{W} \, (\mathbf{O}\mathbf{f}-\mathbf{d}) , \qquad \mathbf{W}=\mathrm{diag}(1/\sigma_k^2) $$
%
% $$ \frac{\partial J}{\partial \mathbf{f}} = \mathbf{O}^T \mathbf{W} (\mathbf{O}\mathbf{f}-\mathbf{d}) , \qquad
%    \frac{\partial^2 J}{\partial \mathbf{f}^2} = \mathbf{O}^T \mathbf{W} \mathbf{O} $$
%
% That is, the finite-element mesh is brought to the data. This is the negative log-likelihood (up to a constant) of the data given f,
% for independent Gaussian errors. There is no normalisation with the area of the domain.
%
% Inputs:
%
%   O     : m x Nnodes observation (node-to-data) operator, as returned by BuildNode2DataMap.m. All rows must correspond to locations
%           inside the mesh, i.e. measurements outside the mesh must have been discarded beforehand (see CreatePointMeasurements.m).  
%   d     : m x 1 measured values
%   dErr  : m x 1 measurement errors (standard deviations), all strictly positive
%   f     : Nnodes x 1 nodal values of the field
%
% Outputs:
%
%   J       : the misfit (a scalar)
%   dJdf    : Nnodes x 1 gradient with respect to the nodal values f
%   d2Jdf2  : Nnodes x Nnodes sparse Hessian, only calculated if requested
%
% Errors are treated as uncorrelated. Spatially correlated errors can be accounted for in an approximate way when preparing the data,
% for example by thinning the data or by using an effective variance, see CreatePointMeasurements.m.
%
% see also: CreatePointMeasurements.m, BuildNode2DataMap.m, Misfit.m, Jq.m, Jqq.m
%
%%

narginchk(4,4)

w=1./dErr(:).^2;
r=O*f(:)-d(:);
wr=w.*r;

J=full(r'*wr)/2;

if nargout>1
    dJdf=full(O'*wr);
end

if nargout>2
    m=numel(w);
    d2Jdf2=O'*spdiags(w,0,m,m)*O;
    d2Jdf2=(d2Jdf2+d2Jdf2.')/2;
end

end
