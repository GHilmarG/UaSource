
%%
function TestPosLinSolution
    
    A=[1 0.1 0.01 0 ; 0.1 1 0 -1 ; 0.1 0 1 0 ; 2 0 -1 1];
   % A=(A+A')/2;
    b=[1 ; 1 ; 1 ; -10000];
    
    I=0;
    iAnew=zeros(size(A,1),1);
    iA=~iAnew;
    while ~all(iAnew==iA) && I<10
        I=I+1;
        fprintf('+++++ Iteration %-i \n ',I)
        iA=iAnew;
        [x,l,iAnew]=SolveConstrainedProblem(A,b,iA);
        x
        l
        if ~all(iAnew==iA)
            fprintf('active set did  change \n')
        else
            fprintf('active set did not change \n')
        end
        %b=b+0.9;
        
    end
    
    [x2,resnorm,residual,exitflag,output,lambda] = lsqnonneg(A,b);
   
    x0=x;
    opts = optimset('Algorithm','active-set','Display','on','MaxIter',100);
    [x3,fval,exitflag] = quadprog(A,-b,[],[],[],[],x*0,[],[],opts);
    
    
    fprintf('%-g \t %-g \t %-g \n',x'*A*x-x'*b,x2'*A*x2-x2'*b,x3'*A*x3-x3'*b)
    fprintf('%-g \t %-g \t %-g \n',norm(A*x-b),norm(A*x2-b),norm(A*x3-b))
    [x x2 x3]
    [A*x-b  A*x2-b A*x3-b]
    
end

function [x,l,iA,c]=SolveConstrainedProblem(A,b,iA)
    
    % on input iA is the current active set
    % on output iA is the new active set
    % if iA on output is not equal iA on input, a feasable (correct) solution has not been found in this iteration
    
    % set up Lagrange matrix
    I=find(iA);
    M=numel(I);
    L=zeros(M,size(A,1)) ; lb=zeros(M,1);
    
    for i=1:M
        L(i,I(i))=1;
    end
    
    Z=zeros(size(L,1),size(L,1));
    AA=[A L' ; L  Z]; bb=[b ; lb];


    sol=AA\bb; 
    x=sol(1:size(A,1)) ; temp=sol(end-size(L,1)+1:end);
    l=zeros(numel(x),1);
    l(I)=temp;
    
    % determine new active set
    % conditions: x_i \ge 0 (pos thickness) , l_i\le 0  , l_i x_i=0 , complementary condition (no summation)
    ix=x>=0; % primal feasibility
    il=l<0;  % A x - b -L' l=0 and l>=0 dual feasibility
    c=ix.*il;  % complementary slackness
    % x=0 iif l<0 (strickt complementary)
    
    % if any component of l is positive then the corresponding element is deleted from the active set 
    % and a new iterate is sought.
    
    iAnew=~ix ;         % these must be included into the new active set
    iAnew(il)=1;        % and these are to be kept in the active est
    iA=iAnew;

end



