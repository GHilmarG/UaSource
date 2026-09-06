function dJ = CalcBruteForceGradient(func,p0,plb,pub,CtrlVar,iRange)

 

%%
%
% p0 : are the model parameters (A,B,C)
%
% func : A function handle
%
%
% [J,G,H]=func(p) 
%
% returns the cost (J), gradient (G) and Hessian (H). (Only cost, J, is used here.) 
%
% The function returns a finite-differences estimate of dJ/dp for selected elements of p in the range iRange
%
% On return dJ is a vector of
%
% [dJ/dA(iRange) dJ/dB(iRange) dJ/dC(iRange)]
%
% assuming one is inverting for A, B and C. 
%
%%

CPUstart=tic;

fprintf(' Calculating gradients using brute-force method...')


J0=func(p0);

%deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize*abs(p0);
deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize*abs(p0)+CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize;


% Testing gradient using brute force method

% deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize*norm(p0);

dJ=p0*0+NaN;
dJtemp=dJ;


parfevalOnAll(gcp('nocreate'), @warning, 0, 'off','MATLAB:decomposition:genericError');
parfevalOnAll(gcp('nocreate'), @warning, 0, 'off','MATLAB:decomposition:SaveNotSupported');

switch  lower(CtrlVar.TestAdjointFiniteDifferenceType)
    
    case {"forward-first-order"}


        parfor k=1:numel(iRange)
            I=iRange(k);
            p1=p0;
            p1(I)=p1(I)+deltaStep(I);


            % what to do if upper or lower bounds are not respected?
            % Here I decide not to try to be too clever, I just test for this and do not calculate if bounds violated
            ind= p1<plb | p1>pub ;

            if any(ind)
                dJtemp(k)=nan;
            else

                J1=func(p1);
                dJtemp(k)=(J1-J0)/deltaStep(I);
            end

        end

        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end

    case {"central-second-order"}


        for k=1:numel(iRange)
            I=iRange(k);
            p1=p0;
            pm1=p0;
            p1(I)=p1(I)+deltaStep(I);
            pm1(I)=pm1(I)-deltaStep(I);


            % what to do if upper or lower bounds are not respected?
            % Here I decide not to try to be too clever, I just test for this and do not calculate if bounds violated
            ind= p1<plb | p1>pub | pm1<plb | pm1 > pub ;

            if any(ind)
                dJtemp(k)=nan;
            else

                J1=func(p1);
                Jm1=func(pm1);
                dJtemp(k)=(J1-Jm1)/deltaStep(I)/2;
            end
        end

        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end


    case {"forward-second-order"}
        
        parfor k=1:numel(iRange)

            I=iRange(k);
            pp1=p0;
            pp2=p0;

            pp1(I)=pp1(I)+deltaStep(I);
            pp2(I)=pp2(I)+2*deltaStep(I);

            % what to do if upper or lower bounds are not respected?
            % Here I decide not to try to be too clever, I just test for this and do not calculate if bounds violated
            ind= pp1<plb | pp1>pub | pp2<plb | pp2 > pub ;

            if any(ind)
                dJtemp(k)=nan;
            else

                Jp1=func(pp1);
                Jp2=func(pp2);
                dJtemp(k)=(-3*J0+4*Jp1-Jp2)/deltaStep(I)/2;
            end

        end
        
        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end

    case "central-fourth-order"

        parfor k=1:numel(iRange)
            I=iRange(k);

            p1=p0;
            pm1=p0;
            p2=p0;
            pm2=p0;
            p1(I)=p1(I)+deltaStep(I);
            p2(I)=p2(I)+2*deltaStep(I);
            pm1(I)=pm1(I)-deltaStep(I);
            pm2(I)=pm2(I)-2*deltaStep(I);

            % what to do if upper or lower bounds are not respected?
            % Here I decide not to try to be too clever, I just test for this and do not calculate if bounds violated
            ind= p1<plb | p1>pub | p2<plb | p2 > pub | pm1 < plb | pm1>pub | pm2 < plb | pm2 > pub; 
            if any(ind)
                dJtemp(k)=nan;
            else

                J1=func(p1);
                J2=func(p2);
                Jm1=func(pm1);
                Jm2=func(pm2);

                dJtemp(k)=(Jm2/12-2*Jm1/3+2*J1/3-J2/12)/deltaStep(I);
            end
            %dJ(I)=(Jm2/12-2*Jm1/3+2*J1/3-J2/12)/deltaStep(I);  % not allowed in parfor
        end

        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end
        
        
    case "complex step differentiation"
        
        fprintf('Error: complex step differentiation not possible because now using a linear solver\n')
        fprintf('that does not work with complex matrices.\n')
        error('CalcBruteForceGradient:CSD','complex step differentiation not possible')
        
        
        parfor k=1:numel(iRange)
            
            I=iRange(k);
            p1=p0;
            p1(I)=p1(I)+1i*deltaStep(I);
            J1=func(p1);
            dJtemp(k)=imag(J1)/deltaStep(I);
            
        end
        
        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end
        
    otherwise
        

        error('which case')
end


CPUend=toc(CPUstart); 
fprintf(' ...done in %f sec.\n',CPUend)
% %%
% calculates d kv/ d b using the complex number method
% using the fact that df/dx=Im(f(x+i dx))/dx
% delta=1e-6*norm(p0);
% dJ=p0*0;
%
% for I=iRange
%
%     p1=p0;
%     p1(I)=p1(I)+1i*delta;
%     J1=func(p1);
%     dJ(I)=imag(J1)/delta;
%
% end

%%
%DkvDb=imag(kv)/db;
%DrhDb=imag(rh)/db;

end


