

%%



CtrlVar=Ua2D_DefaultParameters();
UserVar=[]; 
F=UaFields ; 

F.rhow=1030;

F.h=(0:1000)' ; 
F.rho=F.h*0+920; 
F.S=F.h*0;
F.B=F.h*0-1e10; 

MUA=[];

[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow) ; 


MRP=["0","1","2","3","4"] ;

MRP=["0","l0"] ;
MRP=["1","l1"] ;
MRP=["2","l2"] ;
MRP=["5","l5"] ;
MRP=["3","l3"] ;
MRP=["4","l4"] ;



fig=FindOrCreateFigure("ab") ; clf(fig) ;

hold on

for I=1:numel(MRP)

[F.ab,F.dabdh]=DraftDependentMeltParameterisations(UserVar,CtrlVar,F,MRP(I)) ;

yyaxis left
plot(F.b,F.ab,'b',Displayname="Melt "+MRP(I))
ylabel("basal ablation, $a_b$ (m/yr)",Interpreter="latex")
hold on
yyaxis right 
plot(F.b,F.dabdh,'r',Displayname="Melt derivative "+MRP(I))
ylabel("$d a_b/dh$ (1/yr)",Interpreter="latex")

end

xlabel("lower ice surface, $b$ (m a.s.l.)",Interpreter="latex")

legend(Location="best")

axis padded