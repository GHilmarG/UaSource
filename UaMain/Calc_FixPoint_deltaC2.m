function dIdC=Calc_FixPoint_deltaC2(CtrlVar,UserVar,MUA,F,Meas)

%
%
%
%
%
%  duJdu=(MUA.M*usres)./uErr/Area;
%  dvJdv=(MUA.M*vsres)./vErr/Area;
%  Juv=full(usres'*MUA.M*usres+vsres'*MUA.M*vsres)/2/Area;
%
%  ures=us-Meas.us./uErr
%
% uErr=sqrt(spdiags(Meas.usCov));
%
%   J = ((F.us-Meas.us)./uErr)' * MUA.M *    ((F.us-Meas.us)./uErr) /2/Area
%
%  dJ= M*(u-us)/uErr/Area
%
%  here using u=C*tau
%
%  C tau^m - um


[txzb,tyzb]=CalcNodalStrainRatesAndStresses(CtrlVar,UserVar,MUA,F) ;

tau=sqrt(txzb.*txzb+tyzb.*tyzb) ;



% u=C.*tau.^(m-1).*txzb ;
% v=C.*tau.^(m-1).*tyzb ;

uErr=sqrt(spdiags(Meas.usCov));
vErr=sqrt(spdiags(Meas.vsCov));


Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);

%dIdC=MUA.M*( (F.ub-Meas.us).*txzb./uErr +  (F.vb-Meas.vs).*tyzb./vErr).*(tau.^(m-1)) ./Area ;

% dIdC=MUA.M*( (F.ub-Meas.us).*F.ub./uErr +  (F.vb-Meas.vs).*F.vb./vErr)./F.C./Area ;

dIdC=MUA.M*( (F.ub-Meas.us).*txzb./uErr +  (F.vb-Meas.vs).*tyzb./vErr).*(tau.^(F.m-1)) ./Area ;

% txzb.*tau.^(m-1)  = u/C = F.ub./F.C




% figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdC) ; 

if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
    dIdC=log(10)*F.C.*dIdC;
end

 % dIdC=dIdC/norm(dIdC) ; 
 
end