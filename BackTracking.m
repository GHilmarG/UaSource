function [gmin,fmin,BackTrackInfo,varargout]=BackTracking(slope0,b,fa,fb,Func,CtrlVar,nOut,listInF,listOutF,varargin)

% backtracking to minimize F(x0+gamma DF)
%
% required inputs are slope0, b, fa, fb, F
% if on input slop0 is empty, slope information is not used
% also if on input fa and/or fb are empty they are calculated using a=0 and b=1
% the funcion F must return as a first argument the (scalar) value which is to be minimized
% if nOut is given then further nOut-1 outputs from F are collected in varargout and returned,
% (nOut is then the total number of output variables collected from F)
% the first input to F is the step variable gamma, varargin is then in addition given over to F
%
% F can be called several times. If the output from F are needed as input in a later call to F
% then that can be deal with by defining listInf and litOutF. When F is called
% additional nOut-1 outputs are collected, i.e. the total number of outputs collected from F in nOut
% These additional outputs are then given as additional inputs to F:
%
%
%  [fa,varargout{1:nOut-1}]=F(a,varargin{:}) ;
%   varargin(listInF)=varargout(listOutF)
%
% The general idea behind backtracking is that the search region is shifted
% closer and closer to zero. It is bracketed by [a,b]
% when backtracking:  gamma < b   and  c is the previous b with c > b
% when extrapolating: c=gamma ; b is the previous extrapolation value
%                      b < c=gamma , fb> fc=fgamma
%
%
%
% Example:
%
% 
%%
% 
% 
% 
%   c1=0 ; c2=0.5  ; c3=1.1 ; Func=@(x)  -(x-c1).* (x-c2).*(x-c3) ; slope0=-(c1*c2+c1*c3+c2*c3);  f0=Func(0) ; f1=Func(1) ; 
%
%   [gmin,fmin,BackTrackInfo]=BackTracking(slope0,1,f0,f1,Func);
%
%   xvector=linspace(0,1) ; yvector=Func(xvector);
%   figure
%   plot(xvector,yvector) ; hold on
%   plot(gmin,fmin,'or')
%   plot(BackTrackInfo.Infovector(:,1),BackTrackInfo.Infovector(:,2),'+b')
%
%%
BackTrackInfo.Converged=0;
BackTrackInfo.Direction="   " ;
nFuncEval=0; 
JustPlot=false;

if nargin< 6
    CtrlVar=[];
 
end

if nargin<7
    Fargcollect=0;
    nOut=[] ;
    listInF=[];
    listOutF=[];
else
    Fargcollect=1;
end

if nargin<8
    listInF=[];
    listOutF=[];
end



%%

%% parameters
beta=0.1;
MaxIterations=20;
MaxExtrapolations=100;
ExtrapolationRatio=2.5;
MinXfrac=1e-10*b; % exit if change in x as a fraction of initial interval smaller than this
MaxFuncSame=3; % exit if minimum of func did not change during MaxFuncSame iterations
BacktrackingGammaMin=1e-10;



if isfield(CtrlVar,'BackTrackBeta')
    beta=CtrlVar.BackTrackBeta ;
end
if isfield(CtrlVar,'BackTrackMaxIterations')
    MaxIterations=CtrlVar.BackTrackMaxIterations ;
end
if isfield(CtrlVar,'BackTrackMaxExtrapolations')
    MaxExtrapolations=CtrlVar.BackTrackMaxExtrapolations ;
end
if isfield(CtrlVar,'BackTrackMinXfrac')
    MinXfrac=CtrlVar.BackTrackMinXfrac*b ;
end
if isfield(CtrlVar,'BackTrackMaxFuncSame')
    MaxFuncSame=CtrlVar.BackTrackMaxFuncSame ;
end
if isfield(CtrlVar,'BackTrackExtrapolationRatio')
    ExtrapolationRatio=CtrlVar.BackTrackExtrapolationRatio ;
end

if isfield(CtrlVar,'BacktrackingGammaMin')
    BacktrackingGammaMin=CtrlVar.BacktrackingGammaMin;
end

if ~isfield(CtrlVar,'InfoLevelBackTrack')
    CtrlVar.InfoLevelBackTrack=0;
end

if ~isfield(CtrlVar,'NewtonAcceptRatio')
    CtrlVar.NewtonAcceptRatio=0.5 ;
end

if ~isfield(CtrlVar,'NLtol')
    CtrlVar.NLtol=1e-15;
end

if ~isfield(CtrlVar,'doplots')
    CtrlVar.doplots=0; 
end

if ~isfield(CtrlVar,'BackTrackGuardLower')
    CtrlVar.BackTrackGuardLower=0.45;
end

if ~isfield(CtrlVar,'BackTrackGuardUpper')
    CtrlVar.BackTrackGuardUpper=0.95;
end

if ~isfield(CtrlVar,'LineSearchAllowedToUseExtrapolation')
    CtrlVar.LineSearchAllowedToUseExtrapolation=false;
end

% Backtracking continues even if target has been reached if last reduction in
% ratio is smaller than:
CtrlVar.BackTrackContinueIfLastReductionRatioLessThan=0.5;  

if isempty(slope0) || isnan(slope0)
    NoSlopeInformation=1;
else
    NoSlopeInformation=0;
    if slope0>0
        warning('Backtracking:SlopeAtZeroNotNegative','BackTracking: slope at x=0 must be negative')
        gmin=nan;  fmin=nan; 
        return
    end
end

target=max([CtrlVar.NewtonAcceptRatio*fa CtrlVar.NLtol]);
    
if ~NoSlopeInformation
    targetwithslope=min(fa+beta*slope0*b,CtrlVar.NewtonAcceptRatio*fa);
    target=max(target,targetwithslope);
    
    if target< 0
        %  oops clearly the slope information must be incorrect
        target=max([CtrlVar.NewtonAcceptRatio*fa CtrlVar.NLtol]);
    end
end

xfrac=1e10;

%%
iMinSame=0; iMinSameWhileBacktracking=0;
iarm=0;  BackTrackInfo.iarm=iarm;

a=0;


Infovector=zeros(MaxIterations+3,2)+NaN;



if isempty(fa) || isnan(fa)
    
    if isempty(a) || isnan(a)
        a=0 ;
    end
    
    if Fargcollect
        
        [fa,varargout{1:nOut-1}]=Func(a,varargin{:}) ;
        nFuncEval=nFuncEval+1; 
        
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fa=Func(a);
        nFuncEval=nFuncEval+1; 
    end
    
end

f0=fa; % the value of f at gamma=0

if isempty(fb) || isnan(fb)
    if isempty(b) || isnan(b)
        b=1 ;
    end
    if Fargcollect
        [fb,varargout{1:nOut-1}]=Func(b,varargin{:}) ;
        nFuncEval=nFuncEval+1; 
        
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fb=Func(b);
        nFuncEval=nFuncEval+1; 
    end
end


xStart=b ; fStart=fb; 

Infovector(1:2,1)=[a ; b ] ;  Infovector(1:2,2)=[ fa ; fb ] ; iq=2;  
[fmin,I]=min(Infovector(:,2)) ; gmin=Infovector(I,1) ;

BackTrackInfo.Infovector=Infovector;
BackTrackInfo.nExtrapolationSteps=0;


if fb<target


    if CtrlVar.InfoLevelBackTrack>=10000
        fprintf('B: At start fb<target  (%g<%g). Exiting backtracking \n',fb,target)
    end
    BackTrackInfo.Converged=1;
    BackTrackInfo.nFuncEval=nFuncEval;
    I=isnan(Infovector(:,1)) ; Infovector(I,:)=[];
    BackTrackInfo.Infovector=Infovector;

    if ~isempty(nOut) && nOut> 0
        [fmin,varargout{1:nOut-1}]=Func(gmin,varargin{:}) ;
    end
    % now fmin < ftarget
    % I can now return
    % The only exception is if the user requests information about the
    % backtracking and some pltos

    if ~(CtrlVar.InfoLevelBackTrack>=100 && CtrlVar.doplots==1 )
        return
    else
        JustPlot=true;
    end

end

if ~NoSlopeInformation
    gamma=-b*slope0/2/( (fb-fa)/b-slope0);
else
    %gamma=0.5*a+0.5*b;
    gamma=a+0.9*(b-a);
end

% limit changes and enforce backtracking
if gamma > CtrlVar.BackTrackGuardUpper*b
    gamma=CtrlVar.BackTrackGuardUpper*b ;
elseif gamma < CtrlVar.BackTrackGuardLower*b
    gamma=CtrlVar.BackTrackGuardLower*b;
end


if fmin>target
    iarm=1 ;     BackTrackInfo.iarm=iarm;
    if Fargcollect
        [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{:}) ;
        nFuncEval=nFuncEval+1;
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fgamma=Func(gamma);
        nFuncEval=nFuncEval+1;
    end

    gammaOld=gamma ;

    iq=iq+1 ; Infovector(iq,1)=gamma ;  Infovector(iq,2)=fgamma ;
else
    gammaOld=gamma ;  gamma=b ; fgamma=fb ; 
end

[fmin,I]=min(Infovector(:,2)) ; gmin=Infovector(I,1) ;

c=b; fc=fb ; b=gamma ; fb=fgamma ;

%  I now have values at fa, fb and fc and (possibly) slope0

%if  NoSlopeInformation==1
% target=min([fa,fb])/2;
%else
%    target=fa+beta*slope0*gmin  ; % Armijo criteria
%end

if CtrlVar.InfoLevelBackTrack>=10000
    fprintf('B: step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t f(g)=%-10.5g \t fmin=%-10.5g  \t fmin/ft=%-10.5g \t fmin/f0=%-g \n ',...
        iarm,fa,fb,fc,fgamma,fmin,fmin/target,fmin/f0)
    fprintf('               a=%-10.5g  \t    b=%-10.5g  \t    c=%-10.5g   \t    g=%-10.5g     \t    gmin=%-10.5g \n ',a,b,c,gamma,gmin)
end
%% check extrapolation option

Extrapolation=0;

if CtrlVar.LineSearchAllowedToUseExtrapolation
    while  fb<fa && fc < fb && fmin > target
        Extrapolation=Extrapolation+1;
        
        
        if CtrlVar.InfoLevelBackTrack>=10000
            %    fprintf('Extrapolation step # %-i. fa=%-g \t fb=%-g \t fc=%-g \t fg=%-g \t fmin=%-g \n ',Extrapolation,fa,fb,fc,fgamma,fmin)
            fprintf('E: step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t f(g)=%-10.5g \t fmin=%-10.5g  \t fmin/ft=%-10.5g \t fmin/f0=%-g \n ',...
                Extrapolation-1,fa,fb,fc,fgamma,fmin,fmin/target,fmin/f0)
            fprintf('                a=%-10.5g  \t   b=%-10.5g  \t     c=%-10.5g   \t  g=%-10.5g \t gmin=%-10.5g \n ',a,b,c,gamma,gmin)
            
        end
        
        gamma=ExtrapolationRatio*c ;
        gammaOld=gamma ;
        
        if Fargcollect
            [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
            nFuncEval=nFuncEval+1;
            if ~isempty(listOutF) && ~isempty(listInF)
                [varargin{listInF}]=varargout{listOutF-1} ;
            end
        else
            fgamma=Func(gamma);
            nFuncEval=nFuncEval+1;
        end
        iq=iq+1 ; Infovector(iq,1)=gamma ;  Infovector(iq,2)=fgamma ; [fmin,I]=min(Infovector(:,2)) ; gmin=Infovector(I,1) ;
        
        fbOld=fb; bOld=b;
        fb=fc ; b=c ;
        fc=fgamma ; c=gamma ;
        %
        %         if  ~NoSlopeInformation
        %             target=fa+beta*slope0*gmin  ; % Armijo criteria
        %         end
        %
        
        if Extrapolation > MaxExtrapolations
            if CtrlVar.InfoLevelBackTrack>=10000
                fprintf(' exiting extrapolation step because number of extrapolation steps greater than maximum %-i allowed \n',MaxExtrapolations)
            end
            break
        end
        
    end
end

%%

% if I just came out of extrapolation step, I need to give the final exit results
if CtrlVar.InfoLevelBackTrack>=10000 && Extrapolation>0
    %    fprintf('Extrapolation step # %-i. fa=%-g \t fb=%-g \t fc=%-g \t fg=%-g \t fmin=%-g \n ',Extrapolation,fa,fb,fc,fgamma,fmin)
    fprintf('E: step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t f(g)=%-10.5g \t fmin=%-10.5g  \t f(g)/ft=%-10.5g \t f(g)/f0=%-g \n ',...
        Extrapolation,fa,fb,fc,fgamma,fmin,fgamma/target,fgamma/f0)
    fprintf('                a=%-10.5g  \t b=%-10.5g  \t     c=%-10.5g   \t    g=%-10.5g \t gmin=%-10.5g \n ',a,b,c,gamma,gmin)
    
end

% If extrapolation step was done, then I know that the minimum is to the right of b
% so I shift the 'origin' of the backtrack to b but setting a=b, and then I
% calculate a new value in the middle.
if  Extrapolation>0
    a=bOld ; fa=fbOld; b=(a+c)/2 ;
    NoSlopeInformation=1;
    if Fargcollect
        [fb,varargout{1:nOut-1}]=Func(b,varargin{1:end}) ;
        nFuncEval=nFuncEval+1; 
        
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fb=Func(b);
        nFuncEval=nFuncEval+1; 
    end
    iq=iq+1 ; Infovector(iq,1)=b ;  Infovector(iq,2)=fb ; [fmin,I]=min(Infovector(:,2)) ; gmin=Infovector(I,1) ;
end

fLastReduction=1;
while (fgamma>target || fLastReduction < CtrlVar.BackTrackContinueIfLastReductionRatioLessThan ) && ~JustPlot
    iarm=iarm+1; BackTrackInfo.iarm=iarm;
    
    
    
    if  NoSlopeInformation
        [gamma,pStatus] = parabolamin(a,b,c,fa,fb,fc,CtrlVar.InfoLevelBackTrack);
    else
        [gamma,cStatus]=CubicFit(slope0,fa,fb,fc,b,c,CtrlVar.InfoLevelBackTrack);
        % if cStatus==1
        %     fprintf('Cubic Fit returns status 1 with gamma=%-g \n ',gamma);
        %     [gamma,pStatus] = parabolamin(a,b,c,fa,fb,fc);
        %     if pStatus==1
        %         fprintf('parabolamin returns status 1 with gamma%-g \n ',gamma);
        %     end
        % end
    end
    
    if iarm==2  && Extrapolation>0
        
        % After an extrapolation step
        % I know that the minimum is between a and c with b=(a+c)/2
        
        
        if gamma > b+0.95*(c-b) ; gamma=b+0.95*(c-b) ; elseif gamma < a+0.25*(b-a) ; gamma=a+0.25*(b-a); end
        
        
        if Fargcollect
            [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
            nFuncEval=nFuncEval+1; 
            
            if ~isempty(listOutF) && ~isempty(listInF)
                [varargin{listInF}]=varargout{listOutF-1} ;
            end
        else
            fgamma=Func(gamma);
            nFuncEval=nFuncEval+1; 
        end
        
        
        % Jan 2020: Just changed the limits on gamma so now gamma is somewhere between a and c
        %            But I do not know if gamma is greater or less than b
        
        % I want the minimum to be bracked by a and b
         % b=c ; fb=fc ;
        
        if gamma < b  && fgamma< fb    % min is between a and b
            c=b ; fc=fb ;   % here a < b < c
            b=gamma ; fb=fgamma ;   
        elseif gamma > b && fgamma < fb
            % a=b ; fa=fb ;
            b=c  ; fb=fc ;
            c=gamma ; fc=fgamma ; % now c < b !
        end
        
        %b=gamma ; fb=fgamma ;
        % c=gamma ; fc=fgamma ;
        
    else
        
        % I allow for the possibility that gamma is larger than b and less than c
        % this is the case in the second backstep following extrapolation
        if gamma > b  && gamma < c
            % OK gamma suggested is larger than b
            % so I should try out this value.
            % But make sure the value I try is sufficiently far away from both b and c
            if gamma < (b+CtrlVar.BackTrackGuardUpper*(c-a))
                gamma=b+CtrlVar.BackTrackGuardUpper*(c-a) ;
            elseif gamma > (c-CtrlVar.BackTrackGuardLower*(c-b))
                gamma=c-CtrlVar.BackTrackGuardLower*(c-b);
            end
            
            if Fargcollect
                [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
                nFuncEval=nFuncEval+1;
                
                if ~isempty(listOutF) && ~isempty(listInF)
                    [varargin{listInF}]=varargout{listOutF-1} ;
                end
            else
                fgamma=Func(gamma);
                nFuncEval=nFuncEval+1;
            end
            b=gamma ; fb=fgamma ; % this shifts b to the right
        else
            % general backtracking step
 
            
            if gamma > (a+CtrlVar.BackTrackGuardUpper*(b-a))
                gamma=a+CtrlVar.BackTrackGuardUpper*(b-a) ;
            elseif gamma < (a+CtrlVar.BackTrackGuardLower*(b-a)) 
                gamma=a+CtrlVar.BackTrackGuardLower*(b-a);
            end
            
            
            if Fargcollect
                [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
                nFuncEval=nFuncEval+1; 
                
                if ~isempty(listOutF) && ~isempty(listInF)
                    [varargin{listInF}]=varargout{listOutF-1} ;
                end
            else
                fgamma=Func(gamma);
                nFuncEval=nFuncEval+1; 
            end
            
            c=b; fc=fb ; b=gamma ; fb=fgamma ;  % b always shifted to the left
        end
    end
    %fprintf(' a=%-g \t b=%-g \t c=%-g \t gamma=%-g \n',a,b,c,gamma)
    
    
    fminOld=fmin ; gminOld=gmin;
    iq=iq+1 ; Infovector(iq,1)=gamma ;  Infovector(iq,2)=fgamma ; [fmin,I]=min(Infovector(:,2)) ; gmin=Infovector(I,1) ;
    
    fLastReduction=fmin/fminOld;
    
    if ~isequal(gmin,gminOld) 
        xfrac=abs(gmin-gminOld); 
        MinXfrac=gmin*1e-10 ; % I update this every time a new minimum is found
        
    end
    
    
    if fmin==fminOld
        iMinSame=iMinSame+1;
    else
        iMinSame=0;
    end
    
    if gamma<gammaOld && fmin==fminOld
        iMinSameWhileBacktracking=iMinSameWhileBacktracking+1;
    else
        iMinSameWhileBacktracking=0;
    end
    
    gammaOld=gamma ;
    
    
    if  ~NoSlopeInformation
        target=fa+beta*slope0*b  ; % Armijo criteria
        
        if target< 0
            %  oops clearly the slope information must be incorrect
            target=max([CtrlVar.NewtonAcceptRatio*fa CtrlVar.NLtol]);
        end
    end
    
    %% Print info
    
    if CtrlVar.InfoLevelBackTrack>=10000
        fprintf('B: step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t f(g)=%-10.5g \t fmin=%-10.5g  \t fmin/ft=%-10.5g \t fmin/f0=%-g \n ',...
            iarm,fa,fb,fc,fgamma,fmin,fmin/target,fmin/f0)
        fprintf('                a=%-10.5g  \t    b=%-10.5g  \t    c=%-10.5g  \t   g=%-10.5g \t gmin=%-10.5g \n ',a,b,c,gamma,gmin)
    end
    
    
    %% break criterion
    
    if  iMinSame> MaxFuncSame && fmin < f0
        if CtrlVar.InfoLevelBackTrack>=10000
            fprintf(' exiting backtracking because no further reduction in function over last %-i iterations \n',MaxFuncSame)
        end
        break
    end
    
    if   iarm>5  && iMinSameWhileBacktracking> 2 && fmin < f0
        if CtrlVar.InfoLevelBackTrack>=10000
            fprintf(' exiting backtracking because two subsequent backtracking steps did not result in any further reduction \n')
        end
        break
    end
    
    if xfrac<MinXfrac && fmin < f0
        if CtrlVar.InfoLevelBackTrack>=10000
            fprintf(' exiting backtracking because change in position of minimum %-g less than %-g of interval \n',xfrac,MinXfrac)
        end
        break
    end
    
    if b<BacktrackingGammaMin
        if CtrlVar.InfoLevelBackTrack>=10000
            fprintf(' exiting backtracking because step size (%g) less than minimum allowed step size (%g).\n',b,BacktrackingGammaMin)
        end
        break
    end
    
    
    if iarm>MaxIterations
        if CtrlVar.InfoLevelBackTrack>=10000
            fprintf(' exiting backtracking because number of iteration greater than maximum %-i allowed \n',MaxIterations)
        end
        break
    end
    %%
    
    
end

if fgamma<(target-eps)
    BackTrackInfo.Converged=1;
elseif abs(fmin-fa)/fa > 1e-5 && gamma > BacktrackingGammaMin
    % also consider backtracking to have been success if some reduction was achieved, even though the target reduction may not have been reached
    BackTrackInfo.Converged=1;
else
    BackTrackInfo.Converged=0;  % This is actually value set at start, so not really need, but easy to understand this way
end



%
%     if fgamma<CtrlVar.NewtonAcceptRatio*f0
%         fprintf('B: Fractional reduction (f/f0=%-g) greater than CtrlVar.NewtonAcceptRatio (%-g) \n',fgamma/f0,CtrlVar.NewtonAcceptRatio)
%
%     end
%
%fprintf('B: exit # %-i. fa=%-10.5g \t fb=%-10.5g \t fc=%-10.5g \t fgamma=%-10.5g \t fmin=%-10.5g \t b=%-10.5g \t c=%-10.5g \t gamma=%-10.5g  \n ',Iteration,fa,fb,fc,fgamma,fmin,b,c,gamma)


I=~isnan(Infovector(:,1));  Infovector=Infovector(I,:);
[~,I]=sort(Infovector(:,1)) ; Infovector=Infovector(I,:);

[fgamma,I]=min(Infovector(:,2)) ; gamma=Infovector(I,1) ;



%% Info
if CtrlVar.InfoLevelBackTrack>=100 && CtrlVar.doplots==1

    warning('off','MATLAB:decomposition:SaveNotSupported')
    warning('off','MATLAB:decomposition:genericError')
    parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
    parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:SaveNotSupported');


    if CtrlVar.InfoLevelBackTrack>=1000
        nnn=10 ;

        rTestVector=zeros(nnn,1)+NaN ;
        Upper=1.25*max(Infovector(:,1)) ; Lower=0;
        gammaTestVector=linspace(Lower,Upper,nnn) ;
        dx=min(Infovector(2:end,1)/10) ;
        gammaTestVector=[Lower,dx/1000,dx/50,dx,2*dx,gammaTestVector(2:end)];
        parfor I=1:numel(gammaTestVector)
            gammaTest=gammaTestVector(I);
            rTest=Func(gammaTest);
            gammaTestVector(I)=gammaTest ;
            rTestVector(I)=rTest;
        end
    end

    if isfield(CtrlVar,"BacktracFigName")
        FigName=CtrlVar.BacktracFigName  ;
    else
        FigName="BackTrackingInfo";
    end


    if isfield(CtrlVar,"BacktrackIteration")
        ItText=sprintf("It=%i:  ",CtrlVar.BacktrackIteration);
    else
        ItText="";
    end

    fig=FindOrCreateFigure(FigName) ;  clf(fig) ;
    plot(Infovector(:,1),Infovector(:,2),'or-') ;

    xlabel('$\gamma$',Interpreter='latex') ;
    ylabel('Cost',Interpreter='latex') ;
    
    hold on
    plot(gamma,fgamma,'o',MarkerFaceColor="b",MarkerSize=10)
    
     % add Infovector the the TestVector values

    gammaTestVector=[gammaTestVector(:);Infovector(:,1)];
    rTestVector=[rTestVector(:);Infovector(:,2)];

    [gammaTestVector,Ind]=sort(gammaTestVector) ; rTestVector=rTestVector(Ind) ; 

    if CtrlVar.InfoLevelBackTrack>=1000
       plot(gammaTestVector,rTestVector,'xk-') 
        
    end
    
   
    
    if ~isempty(slope0)
        hold on
        dx=min(Infovector(2:end,1)/10) ;
        plot([0 dx],[f0 f0+slope0*dx],'g','LineWidth',2)
        
    end


    plot(xStart,fStart,"o",MarkerFaceColor="r",MarkerSize=5); 

    legend("backtracking curve values","estimated minimum","cost curve","estimated slope at origin","starting point",Location="best",interpreter="latex")


    title(ItText+sprintf('backtracking/extrapolation steps %-i/%-i',iarm,Extrapolation),Interpreter="latex")

    if isfield(CtrlVar,"time")
        subtitle(sprintf("t=%f   dt=%f",CtrlVar.time,CtrlVar.dt),Interpreter="latex")
    end

    drawnow
    %          prompt = 'Do you want more? Y/N [Y]: ';
    %          str = input(prompt,'s');
    %          if isempty(str)
    %              str = 'Y';
    %          end

end

%%
I=isnan(Infovector(:,1)) ; Infovector(I,:)=[]; 
BackTrackInfo.Infovector=Infovector;
BackTrackInfo.nExtrapolationSteps=Extrapolation;
BackTrackInfo.iarm=iarm;
BackTrackInfo.nFuncEval=nFuncEval;



end
