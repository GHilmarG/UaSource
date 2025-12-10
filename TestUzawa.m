      load TestT A B f g
       
	  [m,n]=size(B);
	  
	  
      y0=g*0;
	  
	  InfoLevel=2 ; IterationMin=3; tol=1e-10;
	  [x,y] = UzawaSymmSolver(A,B,f,g,y0,tol,InfoLevel,IterationMin);
	  
	  
	  