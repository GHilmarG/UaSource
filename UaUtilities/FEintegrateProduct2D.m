



function Int=FEintegrateProduct2D(CtrlVar,MUA,varargin)
%%
%
%  Int=FEintegrateProduct2D(CtrlVar,MUA,varargin)
%
%  calculates the integral of the product of the nodal variable
%  varargin over each of the elements of the FE mesh.
%
%
% Example:
% To calculate the product of h and s over the mesh:
%
%   Int=FEintegrateProduct2D(CtrlVar,MUA,h,s)
%
%
% Testing integration using different quadrature rules of different degrees
% 
%   coordinates=[0 0 ; 0 1 ; 1 1 ; 1 0];
%   connectivity=[1 2  4 ; 4 2 3 ];
%   CtrlVar=Ua2D_DefaultParameters;
%   CtrlVar.TriNodes=3;
%   MUA=CreateMUA(CtrlVar,connectivity,coordinates);
% 
%   CtrlVar.TriNodes=10; MUA=UpdateMUA(CtrlVar,MUA);
%   x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2) ;
%   p=4;
%   f=(1+x+y).^p;
% 
%   I=0 ;
% 
%   F=nan(10,1);
%   D=nan(10,1);
%   for Degree=1:15
% 
% 
%     MUA.QuadratureRuleDegree=Degree;
%     Q=quadtriangle(Degree,'Type','nonproduct','Points','inside','Domain',[0 0 ; 1 0 ; 0 1]) ;
%     MUA.nip=size(Q.Points,1);
%     MUA.niph=size(Q.Points,1);
%     CtrlVar.nip=MUA.nip;
%     CtrlVar.niph=MUA.niph;
%     MUA.points=Q.Points;
%     MUA.weights=Q.Weights;
%     [MUA.Deriv,MUA.DetJ]=CalcMeshDerivatives(CtrlVar,MUA.connectivity,MUA.coordinates,MUA.nip,MUA.points);
% 
%     Int=FEintegrateProduct2D(CtrlVar,MUA,x,y,f,f,f,f,f,f,f,f,f) ;
% 
%     Int=sum(Int);
%     fprintf(' Degree %i \t  nip %i  \t Int=%f \n ',Degree,MUA.nip,Int)
%     I=I+1 ; F(I)=Int; D(I)=Degree; 
% 
%   end
% 
%   figure ; plot(D,F,'-x')
%   xlabel('Degree')
%   ylabel('Fint')
% 
% 
% See Also: FE_inner_product
%
% %%

ndim=2;


Int=zeros(MUA.Nele,1);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    detJw=detJ*MUA.weights(Iint);
    
    fint=1;
    for k=1:numel(varargin)
        
        f=varargin{k};
        fnod=reshape(f(MUA.connectivity,1),MUA.Nele,MUA.nod);
        fint=fnod*fun.*fint;
        
    end
    
    Int=Int+fint.*detJw;
    
end


end