function TestSlope(Func,x,dx,slope)




poolobj = gcp('nocreate'); % If no pool, do not create new one.

if isempty(poolobj)
    parpool(5) ;
else
    poolsize = poolobj.NumWorkers ;
    if poolsize<5
        parpool(5)
    end
end

% spmd
%     F=Func(x+dx*(labindex-3)) ;
%     %labindex
% end

F=zeros(5,1); 
parfor I=-2:2
    F(I+3)=Func(x+dx*I) ;
end

Fm2=F(1);
Fm1=F(2);
F0=F(3) ;
F1=F(4);
F2=F(5) ;

dfdgO1f=(F1-F0)/dx ;                          % first-order forward
dfdgO2f=(-3*F0+4*F1-F2)/dx/2;                 % second-order forward
dfdgO4c=(Fm2/12-2*Fm1/3+2*F1/3-F2/12)/dx;     % fourth-order central

fprintf('         slope       \t\t   first-order forward \t\t  second-order forward    \t fourth-order central \t\t  |slope-dfdfO4)|/|slope|  \n')
fprintf('%15.7g \t\t\t ',slope,dfdgO1f,dfdgO2f,dfdgO4c,norm(slope-dfdgO4c)/norm(slope)) ; fprintf('\n')

%     -logC- after I did the variable change over nodes instead of integration points  !!
%     -2937.486 	    -2937.474 			       -2937.486 			       -2937.486 			    4.074246e-12 nodal-based
%     -2961.843 		-2946.777 			       -2946.789 			       -2946.789 			     0.005082773 int-based
     
% logA  on int, but here A is constant so there is no difference between nodal and int-based values
% -2047.071 			       -2047.065 			       -2047.071 			       -2047.071 			    1.936282e-10 	

% Added sinusoidal pertu
% -2106.471 			       -2166.175 			       -2166.184 			       -2166.176 			       0.0283434 int-based
% -2324.502 			         -2324.5 			       -2324.511 			       -2324.502 			    1.836164e-11 nodal-based

% Another test with -logC- for the whole of Antarctica
% -3.44221e+09 			   -3.247175e+09 			   -3.424744e+09 			   -3.441954e+09 			    7.427023e-05 nodal-based
% -4.746049e+09 			-2.755135e+09 			   -2.984364e+09 			   -3.026078e+09 			       0.3624006 int-based


FindOrCreateFigure('SlopeTest')

plot(dx*[-2 -1 0 1 2],[Fm2,Fm1,F0,F1,F2],'or-')
hold on
plot([0 dx],[F0,F0+slope*dx],'k','LineWidth',2)

end