function [x,y,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g)
    
    
    
    
    if isdiag(B*B')
        
        [nA,nB]=size(A);
        factor=1./(full(mean(abs(diag(A)))));
        A=factor*A ; f=factor*f ;  % this leaves x unaffected but y is scaled
        % factor=1;
        Q=speye(nA,nA)-B'*B;
        
        
        Atilde=Q*A+B'*B;
        btilde=(Q*f+B'*g) ;
        
        
        %dAtilde=factorize(Atilde);
        
        x=Atilde\btilde;
        y=B*(f-A*x);
        
        A=A/factor ; f=f/factor ;  % this leaves x unaffected but y is scaled
        y=y/factor;
        
        tolA=norm(A*x+B'*y-f)/norm(f);
        tolB=norm(B*x-g);
        
        if tolA>1e-6 || tolB>1e-6
            
            fprintf('ABfgPreEliminate: Error seems too large or \t \t \t %g \t %g \n ',norm(A*x+B'*y-f)/norm(f),norm(B*x-g))
            
        end
        
    else
        
        x=NaN;
        y=NaN;
        tolA=NaN;
        tolB=NaN;
        error('ABfgPreEliminate:B','B*B^T not diagonal')
    end
    
    
end


