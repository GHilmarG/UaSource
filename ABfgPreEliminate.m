function [x,y,dAtilde,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g,dAtilde)
    
    
narginchk(6,6)    


    [nA,mA]=size(A);
    [nB,mB]=size(B);
    [nf,mf]=size(f);
   

    if isempty(B) && isempty(g) && ~isempty(A) && ~isempty(f) && mA==nf
        
        % Possibly not needed, but check if this is not just a very simple case of B=g=[]
        
        x=A\f;
        y=NaN;
        
        if nargout>3
            tolB=NaN;
            tolA=norm(A*x-f)/norm(f);
        end
        
    else
        
        BBT=B*B';
        
        if isdiag(BBT)  % the method assumes that B B' is diagonal
            
            
            % It also assumes that B B' is a unity matrix, but if not then simple scaling can be used
            % to ensure that this is the case.
            % To make this a bit more general, I here check if B B' is indeed unity, and
            % if not I do the requried scaling.
            tolerance=eps*1000;
            isBBTunity=all(abs(diag(BBT) - 1) < tolerance) ;
            
            if ~isBBTunity
                [B,g,~,Scale]=ScaleL(CtrlVar,B,g) ;
            else
                Scale=1;
            end
            
            % For numerical reasons a further simple scaling of A is done to bring the
            % sizes of the elements of A in line with those of B.
            factor=1./(full(mean(abs(diag(A)))));
            
            if ~isfinite(factor)  % just in case all elements along the diagonal happen to be equal to zero
                factor=1;
            end
            
            BtB=B'*B;
            
            A=factor*A ; f=factor*f ;  % this leaves x unaffected but y is scaled
            
            
            
            Q=speye(nA,nA)-BtB ; 
            Atilde=Q*A+ BtB ;
            btilde=(Q*f+B'*g) ;

            % https://uk.mathworks.com/help/parallel-computing/benchmarking-a-b.html
            if CtrlVar.Distribute
                if ~isdistributed(Atilde)
                    Atilde=distributed(Atilde);
                end
                if ~isdistributed(btilde)
                    btilde=distributed(btilde);
                end
            end

            %  decomposition is about the same, and as expected this only speeds things up if several solves with the same matrix
            %  are needed.
            %
            tDecomposition=tic;
            if isempty(dAtilde)
                dAtilde=decomposition(Atilde);
            end
            tDecomposition=toc(tDecomposition) ;

            tSolve=tic;
            x=dAtilde\btilde;
            tSolve=toc(tSolve);

            % tSolve=tic;
            % x=Atilde\btilde;
            % tSolve=toc(tSolve);

            [tDecomposition tSolve]

            if isdistributed(x)
                x=gather(x) ;
            end
            % toc

            y=B*(f-A*x);
            
            % Now the solution of the scaled system has been found.
            
            y=y/factor;
            
            if nargout>3
                % check if within tolerances
                
                A=A/factor ;
                f=f/factor ;
                tolA=norm(A*x+B'*y-f)/norm(f);
                tolB=norm(B*x-g);
                
                if tolA>1e-6 || tolB>1e-6
                    
                    fprintf('ABfgPreEliminate: Error seems too large or \t \t \t %g \t %g \n ',norm(A*x+B'*y-f)/norm(f),norm(B*x-g))
                    
                end
            end
            
            y=Scale*y; % and now scale y in case B and g were scaled above.
            
        else
 
 
            error('ABfgPreEliminate:B','B*B^T not diagonal')
        end
        
    end
    
end


