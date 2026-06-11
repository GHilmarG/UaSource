



function Delta=TrustRegionRadiusUpdate(Delta,rho,rhoGood,rhoSmall,rhoTiny,gammai,gammad)


arguments
    Delta  (1,1) {mustBeNumeric}
    rho  (1,1) {mustBeNumeric}
    rhoGood  (1,1) {mustBeNumeric}=0.9
    rhoSmall  (1,1) {mustBeNumeric}=0.5;
    rhoTiny  (1,1) {mustBeNumeric}=0.01;
    gammai  (1,1) {mustBeNumeric}=2;
    gammad (1,1) {mustBeNumeric}=0.5;

end

if rho > rhoGood % very successful, increase Delta
    Delta=gammai*Delta ;
    % elseif rho>etas  ; % not great , but not bad either, keep
    %     Delta=Delta;
elseif rho<rhoTiny
    Delta=Delta/1000; 
elseif rho<rhoSmall  % not successful, decrease delta
    Delta =gammad*Delta;
end



end
