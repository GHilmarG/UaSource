function  [r,rbcs]=ResidualCostFunction(CtrlVar,MUA,L,Rfields,Rbcs,fext0,uvORuvh)
                          


nargoutchk(1,2)
narginchk(7,7)

if isempty(Rbcs) ; Rbcs=0 ; end

rfields=full(abs(Rfields'*Rfields));
rbcs=full(abs(Rbcs'*Rbcs));

f0=full(abs(fext0'*fext0))+1000*eps ; 

r=(rfields+rbcs)/f0;


end

