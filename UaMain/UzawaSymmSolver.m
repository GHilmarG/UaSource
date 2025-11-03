function [x,y] = UzawaSymmSolver(A,B,f,g,y0,CtrlVar)

error('fads')

CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=1; 
[x,y] = AugmentedLagrangianSolver(A,B,f,g,y0,CtrlVar);

return

end

% 
% 
% %[x,y] = UzawaSymmSolver(A,B,f,g,y0,tol,InfoLevel,IterationMin)
% 
% 
% %
% 
% % solves a system on the form [A B' ; B 0] [x; y] =[f;g]
% % using Uzawa iteration, where the inner problem is solved directly
% %
% % A must be symmetrical
% %
% % LDL factorisation done outside of loop, hence cost of each additional iteration fairly small
% %
% 
% % For
% % [A B'] [x]= [f]
% % [B 0 ] [y]  [g]
% % the Uzawa method is:
% %                      A x_{i+1}= f - B' y_i
% %                        y_{i+1}=g y_i + iW B x_{i+1}
% % where iW is `small' but not too small
% 
% %
% % I find that the following modified version of the Uzawa methods converges very quickly:
% % Define T= [A B']
% %           [B iW]
% %
% % Repeat: T [x_{i+1}] = [    f    ]
% %           [y_{i+1]]   [g+ iW y_i]
% %
% % The nice thing is that T has only to be factorized once (LDL if symmetric, LU otherwise) and each additional iteration is cheap
% 
% if CtrlVar.InfoLevelLinSolve>=10
%     fprintf(' Solving a symetrical system with a Uzawa iteration \n' )
% end
% 
% 
% [m,n]=size(B);
% tUzawa=tic;
% 
% % The matrix that I use in the iteration is on the form T=[A B' ; B iW];
% % where iW is `small'.
% % The following seems a good way of setting the values.
% %k=round(log10(max(abs(diag(A))))) ; w=10^(k+CtrlVar.UzawaPower);  iW=speye(m)/w;
% k=round(log10(norm(diag(A))))     ; w=10^(k+CtrlVar.UzawaPower);  iW=speye(m)/w;
% 
% res=1e10 ; res2=1e10 ; 
% T=[A B' ; B iW];
% 
% % LDL
% % P'*S*A*S*P = L*D*L'
% % T=S\(P*L*D*L'*P')/S)
% % S*A*S = P*L*D*L'*P'
% % S=S' , P*P'=1
% % sol=S*P*(L'\(D\(L\(P'*S*fg))));
% 
% 
% tStart=tic;
% %[L,D,P,S]=ldl(T);   % LDL factorisation using MA57
% %O=S*P;  % O'=P'*S'=P'*S
% 
% 
% [L,D,p,S]=ldl(T,'vector');   % LDL factorisation using MA57, MA57 is a multifronta sparse direct solver using AMD ordering
% 
% tElapsed=toc(tStart);
% 
% if CtrlVar.InfoLevelLinSolve >20 ; fprintf(' LDL factorisation in %g sec \n',tElapsed) ; end
% 
% 
% % I measure the residual as res=norm([A B' ; B 0]-[f;g])/norm([x ; y])
% % It is possible that the norm(B*x-g)/norm(y) will become comparable to
% % machine precision and the the iteration should stop, hence the need for res2
% 
% if CtrlVar.InfoLevelLinSolve>1 ;
%     InfoVector=zeros(CtrlVar.UzawaIterationMax,5) ;
% end
% 
% Iteration=0;
% 
% 
% sol=zeros(m+n,1);
%     
% while (res > CtrlVar.LinSolveTol &&  Iteration <= CtrlVar.UzawaIterationMax) || Iteration<  CtrlVar.UzawaIterationMin
%         
%     %tStart=tic;
%     
%     Iteration=Iteration+1;
%     fg=[f ; g + iW*y0];
%     
%     %sol=O*(L'\(D\(L\(O'*fg))));
%     
%     fg=S*fg ; sol(p)=L'\(D\(L\(fg(p)))); sol=S*sol;  % if using the vector format
%     %sol=T\fg;
%     
%     x=sol(1:n) ; y=sol(n+1:end);
%     
%     %[A B'] [x]=[f]
%     %[B 0 ] [y]=[g]
%     
%     res=norm([A B' ; B sparse(m,m)]*sol-[f ; g])/norm([f ;g]);
%     
%     
%     if CtrlVar.InfoLevelLinSolve>1 ;
%         res1=norm(A*x+B'*y-f)/norm(f);
%         ng=norm(g);
%         if ng>0
%             res2=norm(B*x-g)/norm(g);
%         else
%             res2=norm(B*x-g);
%         end
%         
%         InfoVector(Iteration,1)=norm(y-y0)/norm(y); InfoVector(Iteration,2)=res;
%         InfoVector(Iteration,3)=res1; InfoVector(Iteration,4)=res2;
%     end
%     
%     y0=y;
%     
%     
%     
% end
% 
% if CtrlVar.InfoLevelLinSolve>=10
%     disp([' Number of Uzawa2 iterations  ',num2str(Iteration)])
%     tUzawa=toc(tUzawa); disp([' Uzawa  solves in ',num2str(tUzawa),' sec '])
% end
% 
% if Iteration > CtrlVar.UzawaIterationMax
%     fprintf(' Residual %g \n ',res)
%     warning('Uzawa:MaxIterationReached','Uzawa solver exits because maximum number of iterations %g reached \n',CtrlVar.UzawaIterationMax)
% end
% 
% if res > CtrlVar.LinSolveTol 
%     res1=norm(A*x+B'*y-f)/norm(f);
%     res2=norm(B*x-g)/norm(g);
%     fprintf(CtrlVar.fidlog,' Residuals: total  %-g \t first equation %-g \t second equation %-g \n',res,res1,res2);
%     if res2<=10*eps
%         fprintf(CtrlVar.fidlog,' Residuals of second equation comparable to machine precision and iteration therefore stopped \n');
%     else
%     warning('Uzawa:ResidualTooLarge','Uzawa iteration did not fully converge to prescribed tolerance of %-g \n',CtrlVar.LinSolveTol)
%     end
% end
% 
% if CtrlVar.InfoLevelLinSolve>=10
%     fprintf(' ------------ Uzawa symmetrical solver:                                                 --------\n ')
%     fprintf(' Solves a system on the form [A B'' ; B 0] [x;y]=[f;g] , where A=A'', iterativily \n')
%     fprintf(' starting with [A B'' ; B iW]=[f;g]    , where iW=speye(m)/w \n')
%     fprintf(' with w=10^(k+CtrlVar.UzawaPower) , and k=round(log10(norm(diag(A)))) \n')
%     fprintf(' Using k=%g and CtrlVar.UzawaPower=%g \n ',full(k),CtrlVar.UzawaPower)
%     
%     fprintf(' Iteration   y0     norm([A B'' ; B 0]*sol-[f ; g])/norm(sol)     norm(A*x+B''*y-f)/norm([x;y])     norm(B*x-g)/norm(y) \n')
%     for I=1:Iteration
%         fprintf('%10i \t %20.10g \t \t %20.10g \t \t %20.10g \t %20.10g \n',I,InfoVector(I,1),InfoVector(I,2),InfoVector(I,3),InfoVector(I,4))
%     end
%     fprintf(' -----------------------------------------------------------------------------------------------------\n ')
%     iRange=1:Iteration;  figure ; semilogy(iRange,InfoVector(iRange,1),'o-'); xlabel(' Iteration '); ylabel('norm(y)')
%     error('fdsa')
%     
% end
% 
% end
% 
% 
% 
