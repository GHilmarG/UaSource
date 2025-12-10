function q=SchoofFluxFormula(h,A,C,n,m,rho,rhow,g,theta)

narginchk(9,9)

q=rho.* (A*(rho*g).^(n+1).*(1-rho/rhow).^n/(4^n*C^(-1/m))).^(1/(1/m+1)).*theta.^(n/(1/m+1)).*h.^((1/m+n+3)/(1/m+1)) ;

q=real(q);

% units: u rho m = m/yr kg/m^3 m = kg/(m yr) 


end