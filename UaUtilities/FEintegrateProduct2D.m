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
% Testing integration using different number of integration points
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
%   nipv=[1 3 4 6 7 9 12 13 16 19 28 37];
%   F=nipv*0;
%   for nip=nipv
%      CtrlVar.nip=nip ; CtrlVar.niph=nip ;
%      Int=FEintegrateProduct2D(CtrlVar,MUA,x,y,f,f,f,f,f,f,f,f,f) ;
%      Int=sum(Int);
%      fprintf(' nip %i   \t Int=%f \n ',nip,Int)
%      I=I+1 ; F(I)=Int;
%   end
%
%   figure ; plot(nipv,F,'-x')
%   xlabel('nip')
%   ylabel('Fint')
%
%%

ndim=2;

MUA=UpdateMUA(CtrlVar,MUA);

[points,weights]=sample('triangle',MUA.nip,ndim);

Int=zeros(MUA.Nele,1);

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ;
    detJ=MUA.DetJ(:,Iint);
    detJw=detJ*weights(Iint);
    
    fint=1;
    for k=1:numel(varargin)
        
        f=varargin{k};
        fnod=reshape(f(MUA.connectivity,1),MUA.Nele,MUA.nod);
        fint=fnod*fun.*fint;
        
    end
    
    Int=Int+fint.*detJw;
    
end


end