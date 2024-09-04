classdef UaFields


    properties

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
        bmin=[];
        bmax=[];
        bInit=[]

        h=[];
        hInit=[];
        S=[];

        B=[];
        Bmin=[];
        Bmax=[];
        BInit=[];

        AGlen=[];
        AGlenmin=[];
        AGlenmax=[];

        C=[];
        Cmin=[];
        Cmax=[];
        m=[];
        n=[];
        rho=[];
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

        N=[] ;
        aw=[];
        hw=[];
        phi=[] ;
        uw=[];
        vw=[];

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
            if ~isfield(s,'LSF')
                obj.LSF=[];
            end
        end



    end
end