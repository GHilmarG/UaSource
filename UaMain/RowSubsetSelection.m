

function [K,row_idx,flag] = RowSubsetSelection(Aeq, tol)

% ROW_SUBSET_SELECTION  Select p linearly independent rows from L (m x n, rank p <= m)
% If L already has full row rank (p == m), K = L is returned unchanged.
% Otherwise returns K (p x n), full rank, rows are a subset of rows of L.

% check if there is nothing to be done, as Aeq is empty 
if isempty(Aeq)
    K=[];
    flag=0;
    row_idx=[];
    return
end

% if not tolerance is prescribed, set tol to empty 
if nargin < 2
    tol = [];  
end

[m, n] = size(Aeq);


% Column-pivoted QR of L'
[~, R, E] = qr(Aeq', 'vector');

% Estimate rank from diagonal of R
d = abs(diag(R));
if isempty(tol)
    tol = max(m, n) * eps(d(1));  % this is a plausible default tolerance
end
p = full(sum(d > tol));

% fprintf('\t Aeq is %d x %d,  detected rank = %d\n', m, n, p);

% If L has full row rank, no selection needed
if p == m
 %   fprintf('Aeq has full row rank — returning L unchanged.\n');
    K = Aeq;
    flag=0;
    row_idx=[];
    return
end
flag=1; 

% Otherwise select the p most linearly independent rows
row_idx = sort(E(1:p));
K = Aeq(row_idx, :);

% fprintf('\t After row subselection: Aeq is %d x %d,  rank(Aeq) = %d\n', size(K,1), size(K,2), p);
end