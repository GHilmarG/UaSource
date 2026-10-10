




function [UserVar,RunInfo,ub,vb,ud,vd,h]=ExplicitEstimationForUaFields(UserVar,RunInfo,CtrlVar,MUA,F0,Fm1,BCs1,l1,BCs0,l0,dt)
%
% (9 Oct 2026) dt is the new time step, ie F.dt of the field for which the explicit estimate is made. (F0.dt is not
% necessarily the new time step, and CtrlVar.dt is only a copy of F.dt.)
    
    nargoutchk(7,7)
    narginchk(11,11)
    
 
  

    switch CtrlVar.ExplicitEstimationMethod

        case "-no extrapolation-"
           
            h=F0.h ;
            ub=F0.ub;
            vb=F0.vb;
            ud=F0.ud;
            vd=F0.vd;



        case "-dhdt-"



            % [UserVar,dhdt]=dhdtExplicitSUPG(UserVar,CtrlVar,MUA,F0,BCs0);

            [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F0,BCs0); % this now includes rho (2026-04)
            
            h=F0.h+dhdt.*dt ;
            h(h<=CtrlVar.ThickMin)=CtrlVar.ThickMin ;

            % alternative approach: 
            % [UserVar,RunInfo,h]=MassContinuityEquationNewtonRaphsonThicknessContraints(UserVar,RunInfo,CtrlVar,MUA,F0,F0,l0,BCs0) ;
            


            ub=F0.ub+F0.dubdt*dt ;
            vb=F0.vb+F0.dvbdt*dt ;
            ud=F0.ud+F0.duddt*dt ;
            vd=F0.vd+F0.dvddt*dt ;


        case "-Adams-Bashforth-"

            % (9 Oct 2026) The name "-Adams-Bashforth-" is kept for compatibility. The rates in F0 and Fm1 are backward differences
            % over the two previous time steps (see UpdateFtimeDerivatives.m), and the explicit estimate is calculated from these with
            % ExplicitEstimationUsingBackwardDifferences.m, which is second-order accurate also for variable time steps. (Using the
            % backward differences in the Adams-Bashforth formula, as was done previously, is only first-order accurate.) Whether a
            % second-order, linear or no extrapolation can be done is determined from the available data, node by node, within that
            % function. This no longer depends on CtrlVar.CurrentRunStepNumber, and CtrlVar.dtRatio is not used.
            [ub,vb,ud,vd,h]=...
                ExplicitEstimationUsingBackwardDifferences(dt,F0.dtRates,Fm1.dtRates,...
                F0.ub,F0.dubdt,Fm1.dubdt,...
                F0.vb,F0.dvbdt,Fm1.dvbdt,...
                F0.ud,F0.duddt,Fm1.duddt,...
                F0.vd,F0.dvddt,Fm1.dvddt,...
                F0.h,F0.dhdt,Fm1.dhdt);

            % The estimated thickness is not allowed to be below ThickMin. (This is only done for the thickness, not for the velocities.)
            h(h<CtrlVar.ThickMin)=CtrlVar.ThickMin ;

    end


    %%
    %    UaPlots(CtrlVar,MUA,F0,[F0.dubdt F0.dvbdt],FigureTitle="dv/dt")
    %    UaPlots(CtrlVar,MUA,F0,F0.dhdt,FigureTitle="dh/dt")
    %%


    
end
