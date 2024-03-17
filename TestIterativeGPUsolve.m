
%
% As of Matlab 2024a Some GPU operations are much faster. For example, for this test
%
%   DESKTOP-BU2IHIR - GPU NVIDIA Quadro RTX4000
% Matlab 2024a : 1.4sec
% Matlab 2023a : 51.397637 
%
% which is an increase of about factor of 30. This is also explained in the release notes of R2024a
%
%
function TestIterativeGPUsolve

% Create large sparse matrix.
A = gpuArray(gallery("wathen",500,500));

% Create a random right-hand side matrix.
actualSolution = rand(size(A,1),1);
b = A*actualSolution;

% Create preconditioner matrix.
k = 3;
M = tril(triu(A,-k),k);

% Time solver.
maxit = 100;
tol=1e-6; 

tic
[x,flag,relres,iter,resvector]=bicg(A,b,tol,maxit,M) ;
toc

%t=gputimeit(@() bicg(A,b,[],maxNumberOfIterations,M))

figure(10)
semilogy(0:length(resvector)-1,resvector/norm(b),'-o')
yline(tol,'r--');
legend('tri preconditioner','tolerance','Location','East')
xlabel('Iteration number')
ylabel('Relative residual')





end