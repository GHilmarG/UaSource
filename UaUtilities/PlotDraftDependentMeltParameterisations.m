

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


MRP="0" ;

[F.ab,F.dabdh]=DraftDependentMeltParameterisations(UserVar,CtrlVar,F,MRP) ;


fig=FindOrCreateFigure("ab") ; clf(fig) ;

yyaxis left
plot(F.b,F.ab,'b')
ylabel("basal ablation, $a_b$ (m/yr)",Interpreter="latex")
hold on
yyaxis right 
plot(F.b,F.dabdh,'r')
ylabel("$d a_b/dh$ (1/yr)",Interpreter="latex")



xlabel("lower ice surface, $b$ (m a.s.l.)",Interpreter="latex")