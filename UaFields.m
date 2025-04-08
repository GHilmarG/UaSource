classdef UaFields


    properties

        solution="-none-" ; 

        x=[];
        y=[];
        time=[];
        dt=[] ;

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

            obj=s;

            % Make sure the loaded F
            % add in here any new modifications
            if ~isprop(s,'LSF')
                obj.LSF=[];
            end
        end



    end
end