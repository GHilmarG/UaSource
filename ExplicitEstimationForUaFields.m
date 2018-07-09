function F1=ExplicitEstimationForUaFields(CtrlVar,F1,F0,Fm1)

nargoutchk(1,1)
narginchk(4,4)


% The returned F1 does not depend on any of the fields of F1.  I just give this
% as an input so that F1 is not redefined as a structure and to ensure that the
% all other fields are not lost.

%[ub1,vb1,ud1,vd1,h1]=ExplicitEstimation(CtrlVar.dt,dtRatio,CtrlVar.CurrentRunStepNumber,F.ub,F.dubdt,F.dubdtm1,F.vb,F.dvbdt,F.dvbdtm1,F.ud,F.duddt,F.duddtm1,F.vd,F.dvddt,F.dvddtm1,F.h,F.dhdt,F.dhdtm1);


% if the mesh has changed over the last time step, then only derivatives of F0 will have been updated, and not those of Fm1.

if CtrlVar.DebugMode
    filename='Debug_Dumpfile_ExplicitEstimationForUaFields.mat';
   fprintf('ExplicitEstimationForUaFields: Creating dumpfile %s \n',filename) 
end



if CtrlVar.CurrentRunStepNumber>=3
    
    if  (numel(F0.dubdt)~=numel(Fm1.dubdt)) ...
            ||  (numel(F0.dvbdt)~=numel(Fm1.dvbdt)) ...
            ||  (numel(F0.dvbdt)~=numel(Fm1.dvbdt)) ...
            ||  (numel(F0.duddt)~=numel(Fm1.duddt)) ...
            ||  (numel(F0.dvddt)~=numel(Fm1.dvddt))
        
        Fm1.dhdt=F0.dhdt*0;
        Fm1.dubdt=F0.dubdt*0 ; Fm1.dvbdt=F0.dvbdt*0;
        Fm1.duddt=F0.duddt*0 ; Fm1.dvddt=F0.dvddt*0;
    end
    
    CtrlVar.CurrentRunStepNumber=2;
    
end

% improve by checking which fields do need to be updated
[F1.ub,F1.vb,F1.ud,F1.vd,F1.h]=...
    ExplicitEstimation(CtrlVar.dt,CtrlVar.dtRatio,CtrlVar.CurrentRunStepNumber,...
    F0.ub,F0.dubdt,Fm1.dubdt,...
    F0.vb,F0.dvbdt,Fm1.dvbdt,...
    F0.ud,F0.duddt,Fm1.duddt,...
    F0.vd,F0.dvddt,Fm1.dvddt,...
    F0.h,F0.dhdt,Fm1.dhdt);

end
