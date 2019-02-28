
%%

coordinates=[0 0 ; 0 1 ; 1 1 ; 1 0];
connectivity=[1 2  4 ; 4 2 3 ];
CtrlVar=Ua2D_DefaultParameters;
CtrlVar.TriNodes=3;
MUA=CreateMUA(CtrlVar,connectivity,coordinates);

CtrlVar.TriNodes=10; MUA=UpdateMUA(CtrlVar,MUA);
%
%
x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2) ;
p=1;
f=(1+x+y).^p;

I=0 ; 
nipv=[1 3 4 6 7 9 12 13 16 19 28 37];
F=nipv*0;
 for nip=nipv
     CtrlVar.nip=nip ; CtrlVar.niph=nip ;
     Int=FEintegrateProduct2D(CtrlVar,MUA,f,f,f,f,f,f,f,f,f,f,f,f,f);
     %Int=FEintegrate2D(CtrlVar,MUA,f) ;
     Int=sum(Int);
     fprintf(' nip %i   \t Int=%25.15g \n ',nip,Int)
     I=I+1 ; F(I)=Int; 
 end

 figure ; plot(nipv,F,'-x')
 xlabel('nip')
 ylabel('Fint')
