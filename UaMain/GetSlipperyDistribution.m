function [UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F)

narginchk(4,4)
nargoutchk(2,2)

InputFile="DefineSlipperyDistribution.m";
TestIfInputFileInWorkingDirectory(InputFile) ;
NargInputFile=nargin(InputFile);
N=nargout('DefineSlipperyDistribution');



SlidingLaw=SlidingLawCode(CtrlVar.SlidingLaw);   % (26 Oct 2026) names and aliases, see SlidingLawCode.m

if SlidingLaw==2 && N<4   % Budd
    error("GetSlipperyDistribution:nargout","When using Budd sliding law, DefineSlipperyDistribution.m must return 4 parameters [UserVar,C,m,q] ")
end


if any(SlidingLaw==[3 4 5 6]) && N<5   % Coulomb, Tsai, Cornford, Umbi (need muk)
    fprintf("\n \n When using the sliding law %s, DefineSlipperyDistribution.m must return 5 parameters [UserVar,C,m,q,muk]. \n",CtrlVar.SlidingLaw)
    fprintf(" Note: The sliding law only depends on the parameters C, m, muk. (so you can set, for example, q=NaN.)  \n")
    error("GetSlipperyDistribution:nargout","Incorrect number of output parameters in DefineSlipperiness ")
end


switch N

    case 3

        if NargInputFile>4

            [UserVar,F.C,F.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);

        else

            [UserVar,F.C,F.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,F);
        end

    case 4

        if NargInputFile>4
            [UserVar,F.C,F.m,F.q]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);
        else
            [UserVar,F.C,F.m,F.q]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,F);
        end

    case 5

        if NargInputFile>4

            [UserVar,F.C,F.m,F.q,F.muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);

        else

            [UserVar,F.C,F.m,F.q,F.muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,F);

        end

    case 6

        if NargInputFile>4

            [UserVar,F.C,F.m,F.q,F.muk,F.V0]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF);

        else

            [UserVar,F.C,F.m,F.q,F.muk,F.V0]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,F);

        end

    otherwise

        error('Ua:GetSlipperyDistribution','DefineSlipperyDistribution must return between 3 and 5 arguments')

end

F.Cmax=CtrlVar.Cmax;
F.Cmin=CtrlVar.Cmin;

[F.C,F.m,F.q,F.muk,F.V0]=TestSlipperinessInputValues(CtrlVar,MUA,F.C,F.m,F.q,F.muk,F.V0);


if CtrlVar.IncludeMelangeModelPhysics

    [UserVar,F]=GetSeaIceParameters(UserVar,CtrlVar,MUA,F) ;


end

end




