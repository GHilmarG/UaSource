




function [UserVar,RunInfo,ub,vb,ud,vd,h]=ExplicitEstimationForUaFields(UserVar,RunInfo,CtrlVar,MUA,F0,Fm1,BCs1,l1,BCs0,l0)
    
    nargoutchk(7,7)
    narginchk(10,10)
    
 
  

    switch CtrlVar.ExplicitEstimationMethod

        case "-no extrapolation-"
           
            h=F0.h ;
            ub=F0.ub;
            vb=F0.vb;
            ud=F0.ud;
            vd=F0.vd;



        case "-dhdt-"



            [UserVar,dhdt]=dhdtExplicitSUPG(UserVar,CtrlVar,MUA,F0,BCs0);
            h=F0.h+dhdt.*CtrlVar.dt ;
            h(h<=CtrlVar.ThickMin)=CtrlVar.ThickMin ;

            % alternative approach: 
            % [UserVar,RunInfo,h]=MassContinuityEquationNewtonRaphsonThicknessContraints(UserVar,RunInfo,CtrlVar,MUA,F0,F0,l0,BCs0) ;
            


            ub=F0.ub+F0.dubdt*CtrlVar.dt ;
            vb=F0.vb+F0.dvbdt*CtrlVar.dt ;
            ud=F0.ud+F0.duddt*CtrlVar.dt ;
            vd=F0.vd+F0.dvddt*CtrlVar.dt ;


        case "-Adams-Bashforth-"

            
            
            Itime=CtrlVar.CurrentRunStepNumber;
            
            if CtrlVar.CurrentRunStepNumber>=3
                
                if     (numel(F0.ub)~=numel(F0.dubdt)) ...
                        || (numel(F0.vb)~=numel(F0.dvbdt)) ...
                        || (numel(F0.ud)~=numel(F0.duddt)) ...
                        || (numel(F0.vd)~=numel(F0.dvddt))
                    
                    Itime=1;
                    
                elseif      (numel(F0.dubdt)~=numel(Fm1.dubdt)) ...
                        ||  (numel(F0.dvbdt)~=numel(Fm1.dvbdt)) ...
                        ||  (numel(F0.dvbdt)~=numel(Fm1.dvbdt)) ...
                        ||  (numel(F0.duddt)~=numel(Fm1.duddt)) ...
                        ||  (numel(F0.dvddt)~=numel(Fm1.dvddt))
                    Itime=2 ;
                end
                
            end
            
            
            
         
            % improve by checking which fields do need to be updated
            [ub,vb,ud,vd,h]=...
                ExplicitEstimation(CtrlVar.dt,CtrlVar.dtRatio,Itime,...
                F0.ub,F0.dubdt,Fm1.dubdt,...
                F0.vb,F0.dvbdt,Fm1.dvbdt,...
                F0.ud,F0.duddt,Fm1.duddt,...
                F0.vd,F0.dvddt,Fm1.dvddt,...
                F0.h,F0.dhdt,Fm1.dhdt);

    end


    %%
    %    UaPlots(CtrlVar,MUA,F0,[F0.dubdt F0.dvbdt],FigureTitle="dv/dt")
    %    UaPlots(CtrlVar,MUA,F0,F0.dhdt,FigureTitle="dh/dt")
    %%


    
end
