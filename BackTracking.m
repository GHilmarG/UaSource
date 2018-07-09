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
%  [fa,varargout{1:nOut-1}]=F(a,varargin{:}) ;
%   varargin(listInF)=varargout(listOutF)
%
% The general idea behind backtracking is that the search region is shifted
% closer and closer to zero. It is bracketed by [a,b]
% when backtracking:  gamma < b   and  c is the previous b with c > b
% when extrapolating: c=gamma ; b is the previous extrapolation value
%                      b < c=gamma , fb> fc=fgamma



%%
BackTrackInfo.converged=0;
%% First check if the input value can already be accepted

target=max([CtrlVar.NewtonAcceptRatio*fa CtrlVar.NLtol]);
xfrac=1e10;

if nargin<7
    Fargcollect=0;
else
    Fargcollect=1;
end


if fb<target
    
    gamma=b; fgamma=fb; Iteration=0; fc=NaN; fmin=fb ; gmin=b ; c=NaN ;
    InfoVector(1:2,1)=[0 ; b ] ;  InfoVector(1:2,2)=[ fa ; fb ] ;
    BackTrackInfo.InfoVector=InfoVector;
    BackTrackInfo.nExtrapolatinSteps=0;
    BackTrackInfo.nBackTrackSteps=Iteration;
    
    
    if Fargcollect
        
        [fb,varargout{1:nOut-1}]=Func(b,varargin{:}) ;
        
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fb=Func(b);
    end
    
    %[fb,varargout{1:nOut-1}]=F(b,varargin{:}) ;
    
    
    if CtrlVar.InfoLevelBackTrack>=2
        fprintf('B: At start fb<target  (%g<%g). Exiting backtracking \n',fb,target)
    end
    BackTrackInfo.converged=1;
    return
end

%%




%%


%% parameters
beta=0.1;
MaxIterations=50;
MaxExtrapolations=100;
ExtrapolationRatio=2.5;
MinXfrac=1e-10*b; % exit if change in x as a fraction of initial interval smaller than this
MaxFuncSame=3; % exit if minimum of func did not change during MaxFuncSame iterations
BacktrackingGammaMin=1e-20;
InfoLevelBackTrack=1;

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

if isfield(CtrlVar,'InfoLevelBackTrack')
 InfoLevelBackTrack=CtrlVar.InfoLevelBackTrack;
end
%%
iMinSame=0; iMinSameWhileBacktracking=0;
a=0;



InfoVector=zeros(MaxIterations+3,2)+NaN;

if isempty(slope0)
    NoSlopeInformation=1;
else
    NoSlopeInformation=0;
    if slope0>0
        error('BackTracking: slope at x=0 must be negative')
    end
end


if isempty(fa)
    a=0 ;
    
    if Fargcollect
        
        [fa,varargout{1:nOut-1}]=Func(a,varargin{:}) ;
        
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fa=Func(a);
    end
    
end



f0=fa; % the value of f at gamma=0

if isempty(fb)
    b=1 ;
    
    if Fargcollect
        [fb,varargout{1:nOut-1}]=Func(b,varargin{:}) ;
        
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fb=Func(b);
    end
    
    
end



InfoVector(1:2,1)=[a ; b ] ;  InfoVector(1:2,2)=[ fa ; fb ] ; iq=2;  [fmin,I]=min(InfoVector(:,2)) ; gmin=InfoVector(I,1) ;


if ~NoSlopeInformation
    gamma=-b*slope0/2/( (fb-fa)/b-slope0);
else
    %gamma=0.5*a+0.5*b;
    gamma=a+0.9*(b-a);
end

% limit changes and enforce backtracking
if gamma > 0.9*b ; gamma=0.9*b ; elseif gamma < 0.25*b ; gamma=0.25*b; end

Iteration=1 ;
if Fargcollect
    [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{:}) ;
    if ~isempty(listOutF) && ~isempty(listInF)
        [varargin{listInF}]=varargout{listOutF-1} ;
    end
else
    fgamma=Func(gamma);
end

gammaOld=gamma ;

iq=iq+1 ; InfoVector(iq,1)=gamma ;  InfoVector(iq,2)=fgamma ; [fmin,I]=min(InfoVector(:,2)) ; gmin=InfoVector(I,1) ;

c=b; fc=fb ; b=gamma ; fb=fgamma ;

%  I now have values at fa, fb and fc and (possibly) slope0

%if  NoSlopeInformation==1
target=min([fa,fb])/2;
%else
%    target=fa+beta*slope0*gmin  ; % Armijo criteria
%end

if InfoLevelBackTrack>=2
    fprintf('B: step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t f(g)=%-10.5g \t fmin=%-10.5g  \t fmin/ft=%-10.5g \t fmin/f0=%-g \n ',...
        Iteration,fa,fb,fc,fgamma,fmin,fmin/target,fmin/f0)
    fprintf('               a=%-10.5g  \t    b=%-10.5g  \t    c=%-10.5g   \t    g=%-10.5g     \t    gmin=%-10.5g \n ',a,b,c,gamma,gmin)
end
%% check extrapolation option

Extrapolation=0;
while  fb<fa && fc < fb && fmin > target
    Extrapolation=Extrapolation+1;
    
    
    if InfoLevelBackTrack>=2
        %    fprintf('Extrapolation step # %-i. fa=%-g \t fb=%-g \t fc=%-g \t fg=%-g \t fmin=%-g \n ',Extrapolation,fa,fb,fc,fgamma,fmin)
        fprintf('E: step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t f(g)=%-10.5g \t fmin=%-10.5g  \t fmin/ft=%-10.5g \t fmin/f0=%-g \n ',...
            Extrapolation-1,fa,fb,fc,fgamma,fmin,fmin/target,fmin/f0)
        fprintf('                a=%-10.5g  \t   b=%-10.5g  \t     c=%-10.5g   \t  g=%-10.5g \t gmin=%-10.5g \n ',a,b,c,gamma,gmin)
        
    end
    
    gamma=ExtrapolationRatio*c ;
    gammaOld=gamma ;
    
    if Fargcollect
        [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fgamma=Func(gamma);
    end
    iq=iq+1 ; InfoVector(iq,1)=gamma ;  InfoVector(iq,2)=fgamma ; [fmin,I]=min(InfoVector(:,2)) ; gmin=InfoVector(I,1) ;
    
    fbOld=fb; bOld=b;
    fb=fc ; b=c ;
    fc=fgamma ; c=gamma ;
    %
    %         if  ~NoSlopeInformation
    %             target=fa+beta*slope0*gmin  ; % Armijo criteria
    %         end
    %
    
    if Extrapolation > MaxExtrapolations
        if InfoLevelBackTrack>=2
            fprintf(' exiting extrapolation step because number of extrapolation steps greater than maximum %-i allowed \n',MaxExtrapolations)
        end
        break
    end
    
end


%%

% if I just came out of extrapolation step, I need to give the final exit results
if InfoLevelBackTrack>=2 && Extrapolation>0
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
        
        if ~isempty(listOutF) && ~isempty(listInF)
            [varargin{listInF}]=varargout{listOutF-1} ;
        end
    else
        fb=Func(b);
    end
    iq=iq+1 ; InfoVector(iq,1)=b ;  InfoVector(iq,2)=fb ; [fmin,I]=min(InfoVector(:,2)) ; gmin=InfoVector(I,1) ;
end

fLastReduction=1;
while fgamma>target || fLastReduction < 0.5
    Iteration=Iteration+1;
    
    
    
    if  NoSlopeInformation
        [gamma,pStatus] = parabolamin(a,b,c,fa,fb,fc,InfoLevelBackTrack);
    else
        [gamma,cStatus]=CubicFit(slope0,fa,fb,fc,b,c,InfoLevelBackTrack);
        if cStatus==1
            fprintf(CtrlVar.fidlog,'Cubic Fit returns status 1 with gamma=%-g \n ',gamma);
            [gamma,pStatus] = parabolamin(a,b,c,fa,fb,fc);
            if pStatus==1
                fprintf(CtrlVar.fidlog,'parabolamin returns status 1 with gamma%-g \n ',gamma);
            end
        end
    end
    
    if Iteration==2  && Extrapolation>0
        
        % After an extrapolation step
        % I know that the minimum is between a and c with b=(a+c)/2
        
        if gamma > b+0.95*(c-b) ; gamma=b+0.95*(c-b) ; elseif gamma < b+0.75*(c-b) ; gamma=b+0.75*(c-b); end
        %fprintf(' a=%-g \t b=%-g \t c=%-g \t gamma=%-g \n',a,b,c,gamma)
        
        if Fargcollect
            [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
            
            if ~isempty(listOutF) && ~isempty(listInF)
                [varargin{listInF}]=varargout{listOutF-1} ;
            end
        else
            fgamma=Func(gamma);
        end
        
        
        %b=gamma ; fb=fgamma ;
        c=gamma ; fc=fgamma ;
    else
        
        % I allow for the possibility that gamma is larger than b and less than c
        % this is the case in the second backstep following extrapolation
        if gamma >=b+0.5*(c-b) && gamma <=b+0.95*(c-b)
            if Fargcollect
                [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
                
                if ~isempty(listOutF) && ~isempty(listInF)
                    [varargin{listInF}]=varargout{listOutF-1} ;
                end
            else
                fgamma=Func(gamma);
            end
            b=gamma ; fb=fgamma ; % this shifts b to the right
        else
            % general backtracking step
            %  a+0.25 (b-a)  < gamma < a+ 0.95 (b-a)
            if gamma > (a+0.95*(b-a)) ; gamma=a+0.95*(b-a) ; elseif gamma < (a+0.25*(b-a)) ; gamma=a+0.25*(b-a); end
            if Fargcollect
                [fgamma,varargout{1:nOut-1}]=Func(gamma,varargin{1:end}) ;
                
                if ~isempty(listOutF) && ~isempty(listInF)
                    [varargin{listInF}]=varargout{listOutF-1} ;
                end
            else
                fgamma=Func(gamma);
            end
            
            c=b; fc=fb ; b=gamma ; fb=fgamma ;  % b always shifted to the left
        end
    end
    %fprintf(' a=%-g \t b=%-g \t c=%-g \t gamma=%-g \n',a,b,c,gamma)
    
    
    fminOld=fmin ; gminOld=gmin;
    iq=iq+1 ; InfoVector(iq,1)=gamma ;  InfoVector(iq,2)=fgamma ; [fmin,I]=min(InfoVector(:,2)) ; gmin=InfoVector(I,1) ;
    
    fLastReduction=fmin/fminOld;
    if ~isequal(gmin,gminOld) ; xfrac=abs(gmin-gminOld); end
    
    
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
    end
    
    if  iMinSame> MaxFuncSame && fmin < f0
        if InfoLevelBackTrack>=2
            fprintf(' exiting backtracking because no further reduction in function over last %-i iterations \n',MaxFuncSame)
        end
        break
    end
    
    if  iMinSameWhileBacktracking> 4 && fmin < f0
        if InfoLevelBackTrack>=2
            fprintf(' exiting backtracking because two subsequent backtracking steps did not result in any further reduction \n')
        end
        break
    end
    
    if xfrac<MinXfrac
        if InfoLevelBackTrack>=2
            fprintf(' exiting backtracking because change in position of minimum %-g less than %-g of interval \n',xfrac,MinXfrac)
        end
        break
    end
    
    if b<BacktrackingGammaMin
        if CtrlVar.InfoLevelBackTrack>=2
            fprintf(' exiting backtracking because step size (%g) less than minimum allowed step size (%g).\n',b,BacktrackingGammaMin)
        end
        break
    end
    
    
    if Iteration>MaxIterations
        if InfoLevelBackTrack>=2
            fprintf(' exiting backtracking because number of iteration greater than maximum %-i allowed \n',MaxIterations)
        end
        break
    end
    
    if InfoLevelBackTrack>=2
        fprintf('B: step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t f(g)=%-10.5g \t fmin=%-10.5g  \t fmin/ft=%-10.5g \t fmin/f0=%-g \n ',...
            Iteration,fa,fb,fc,fgamma,fmin,fmin/target,fmin/f0)
        fprintf('                a=%-10.5g  \t    b=%-10.5g  \t    c=%-10.5g  \t   g=%-10.5g \t gmin=%-10.5g \n ',a,b,c,gamma,gmin)
    end
    
end

if fmin<fa
    BackTrackInfo.converged=1;
end



%
%     if fgamma<CtrlVar.NewtonAcceptRatio*f0
%         fprintf('B: Fractional reduction (f/f0=%-g) greater than CtrlVar.NewtonAcceptRatio (%-g) \n',fgamma/f0,CtrlVar.NewtonAcceptRatio)
%
%     end
%
%fprintf('B: exit # %-i. fa=%-10.5g \t fb=%-10.5g \t fc=%-10.5g \t fgamma=%-10.5g \t fmin=%-10.5g \t b=%-10.5g \t c=%-10.5g \t gamma=%-10.5g  \n ',Iteration,fa,fb,fc,fgamma,fmin,b,c,gamma)


I=~isnan(InfoVector(:,1));  InfoVector=InfoVector(I,:);
[~,I]=sort(InfoVector(:,1)) ; InfoVector=InfoVector(I,:);

[fgamma,I]=min(InfoVector(:,2)) ; gamma=InfoVector(I,1) ;

if InfoLevelBackTrack>=100 && CtrlVar.doplots==1
    
    figure
    plot(InfoVector(:,1),InfoVector(:,2),'or-') ; xlabel('gamma') ; ylabel('Cost') ;
    title(sprintf('backtracking/extrapolation steps %-i/%-i',Iteration,Extrapolation))
    hold on
    plot(gamma,fgamma,'xg')
    
    if ~NoSlopeInformation
        hold on
        dx=min(InfoVector(2:end,1)/10) ;
        plot([0 dx],[f0 f0+slope0*dx],'g')
        hold off
    end
    
    %          prompt = 'Do you want more? Y/N [Y]: ';
    %          str = input(prompt,'s');
    %          if isempty(str)
    %              str = 'Y';
    %          end
end


BackTrackInfo.InfoVector=InfoVector;
BackTrackInfo.nExtrapolatinSteps=Extrapolation;
BackTrackInfo.nBackTrackSteps=Iteration;



end
