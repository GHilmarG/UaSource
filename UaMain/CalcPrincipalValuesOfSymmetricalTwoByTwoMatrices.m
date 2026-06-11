







function [pV,eV1,eV2]=CalcPrincipalValuesOfSymmetricalTwoByTwoMatrices(txx,txy,tyy,options)

%%
%
%   T=[txx  txy]
%     [txy  tyy]
%
%
%
% Outputs:
%
%   pV   : principal strain rates, pV(:,1) is the first, and pV(:,2) the second.
%
%%

arguments

    txx      (:,1) double
    txy      (:,1) double
    tyy     (:,1) double
    options.eigenvalues logical = true
    options.eigenvectors logical = false

end

N=numel(txx);
pV=zeros(N,2);

eV1=[];
eV2=[];





if options.eigenvalues && ~options.eigenvectors

    % If only the eigenvalues are needed, presumably quickest to calculate this directly
    % The way I calculate this, pV(:,2) will always be the larger eigenvalue of the two.


    SQ=sqrt( 4*txy.*txy + (txx-tyy).^2) ;
    pV(:,1)=( (txx+tyy) - SQ )/2;
    pV(:,2)=( (txx+tyy) + SQ )/2;


else

    eV1=zeros(N,2);
    eV2=zeros(N,2);

    for k=1:N

        T=[txx(k) txy(k) ; txy(k) tyy(k)];

        [V,D]=eig(T);

        pV(k,1)=D(1,1);
        pV(k,2)=D(2,2);
        eV1(k,1)=V(1,1) ;   % x value of first eigen vector
        eV1(k,2)=V(2,1) ;   % y value of first eigen vector
        eV2(k,1)=V(1,2) ;   % x value of second eigen vector
        eV2(k,2)=V(2,2) ;   % y value of second eigen vector


    end

end