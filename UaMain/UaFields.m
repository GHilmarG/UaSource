classdef UaFields


    properties

        solution="-none-" ;

        x=[];
        y=[];
        time=[];
        dt=[] ;
        dtRates=NaN ;  % (9 Oct 2026) the time step over which the rates (dhdt, dubdt, dvbdt, duddt, dvddt) were calculated as backward differences

        hTimeDiscretisationErrorEstimate=[] ;  % (10 Oct 2026) estimate of the local time-discretisation error of h in the last time step (m),
        % see TimeDiscretisationErrorEstimate.m

        hTimeDiscretisationErrorAccumulated=[] ;  % (10 Oct 2026) accumulated estimate (structure with fields Signed, Abs, ValidTime, TotalTime),
        % see TimeDiscretisationErrorEstimate.m

        ah=[] ;  % (15 Oct 2026) equivalent local mass-balance rate of the thickness constraints (distance/time), positive where ice is added,
        % see ThicknessConstraintReactions.m

        lGL=[] ;  % (2 Nov 2026) nodal transition length of the grounding-line mask (m), or [] for the scalar CtrlVar.GroundingLineTransition.lmin,
        % set at the start of each time step, see GroundingLineTransitionLength.m

        hSuppressedThickening=[] ;  % (31 Oct 2026) net thickening suppressed by the thickness constraint at each constrained node (distance),
        % for the cumulative release criterion, see UpdateSuppressedThickening.m

        xint=[];
        yint=[];

        ub=[];
        vb=[];

        ud=[];
        vd=[];

        uo=[];
        vo=[];

        ua=[];
        va=[];



        s=[] ;
        sInit=[];

        b=[];
        b0=[];
        bmin=[];
        bmax=[];
        bInit=[]

        h=[];
        h0=[];
        hInit=[];
        S=[];

        B=[];
        Bmin=[];
        Bmax=[];
        BInit=[];

        AGlen=[];
        AGlenmin=[];
        AGlenmax=[];
        AGlen0=[] ; % undamaged A, used in phase field fracture

        E=[];

        C=[];
        Cmin=[];
        Cmax=[];

        m=[];
        n=[];
        rho=[];
        rho0=[]; % undamaged rho, used in phase field fracture
        rhow=[];

        q=[];
        muk=[];
        V0=[] ; % This is a parameter in Joughin's sliding law, rCW-V0

        Co=[];
        mo=[]
        Ca=[];
        ma=[];

        as=[];
        ab=[];
        dasdh=[];
        dabdh=[];


        dhdt=[] ;
        dsdt=[] ;
        dbdt=[] ;

        dubdt=[];
        dvbdt=[];

        duddt=[];
        dvddt=[];



        g=[];
        alpha=0;



        GF=[];
        GFInit=[];

        LSF=[] % Level Set Field
        LSFMask=[];
        LSFnodes=[];
        c=[] ; % calving rate
        LSFqx=[] ;
        LSFqy=[] ;

        % subglacier water
        N=[] ;
        aw=[];
        hw=[];
        phi=[] ;  % also used for phase field fracture
        uw=[];
        vw=[];


        Psi=[] ; % strain-rate energy density function

        D=[] ; % Damage (SSD)
        aD=[] ; % Damage accumulation

        txx=[];
        txy=[];
        tyy=[];


    end



    methods (Static)

        function obj = loadobj(s)

            % (2 Nov 2026) If the object was saved with a different set of properties (eg with a property that has since been removed,
            % such as kHElement), MATLAB passes it as a structure. A UaFields object is then created and the fields that are still
            % properties are copied; fields that no longer exist are dropped. (Returning the structure itself gave the warning
            % MATLAB:class:mustReturnObject.)
            if isstruct(s)
                obj=UaFields ;
                fn=intersect(fieldnames(s),properties(obj)) ;
                for k=1:numel(fn)
                    obj.(fn{k})=s.(fn{k}) ;
                end
            else
                obj=s;
            end

            % Make sure the loaded F
            % add in here any new modifications
            if ~isprop(obj,'LSF')
                obj.LSF=[];
            end
        end



    end
end