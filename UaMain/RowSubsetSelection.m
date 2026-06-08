

function [K,row_idx,flag] = RowSubsetSelection(L, tol)

% ROW_SUBSET_SELECTION  Select p linearly independent rows from L (m x n, rank p <= m)
% If L already has full row rank (p == m), K = L is returned unchanged.
% Otherwise returns K (p x n), full rank, rows are a subset of rows of L.

if nargin < 2
    tol = [];   % auto-detect
end

[m, n] = size(L);

% Column-pivoted QR of L'
[~, R, E] = qr(L', 'vector');

% Estimate rank from diagonal of R
d = abs(diag(R));
if isempty(tol)
    tol = max(m, n) * eps(d(1));
end
p = full(sum(d > tol));

fprintf('L is %d x %d,  detected rank = %d\n', m, n, p);

% If L has full row rank, no selection needed
if p == m
    fprintf('L has full row rank — returning L unchanged.\n');
    K = L;
    flag=0;
    return
end
flag=1; 

% Otherwise select the p most linearly independent rows
row_idx = sort(E(1:p));
K = L(row_idx, :);

fprintf('K is %d x %d,  rank(K) = %d\n', size(K,1), size(K,2), p);
end